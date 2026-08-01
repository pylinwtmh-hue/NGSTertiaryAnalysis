#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=========================================================
WGS/WES Germline Tertiary Analysis - build_clingen_erepo_lookup.py
=========================================================
Author : Po-Yu Lin (林伯昱)  <p88124019@gs.ncku.edu.tw>
Licensed under the GNU General Public License v3.0
=========================================================
把 ClinGen Evidence Repository（ERepo）的全量下載檔轉成 pipeline 用的查表。

ERepo 是什麼：ClinGen 各 Variant Curation Expert Panel（VCEP）對變異做的**專家判讀**，
不只給結論（Pathogenic / Likely Pathogenic / VUS / …），還記錄「實際套用了哪些 ACMG
criteria」（例如 PVS1, PM2_Supporting, PP3）。這正好可以拿來跟我們自動 ACMG 的結果對照。

⚠️ 定位：**只作對照（QC），不進 ACMG 計分**。
   理由：ClinGen SVI 2018 已建議不要用 PP5/BP6（拿他人判讀當證據），否則會變成循環論證。
   我們的 pipeline 也刻意沒有實作 PP5/BP6，這裡維持同一原則。
   註：VCEP 判讀本來就會回存 ClinVar，並顯示為 review status「reviewed by expert panel」
   （3 星）；我們的 CLINVAR_STARS 已能看到「有沒有專家判讀」，ERepo 補的是「用了哪些
   criteria」這一層細節。

join key：**ClinVar Variation ID**（ERepo 以 ClinVar Variation ID 或 ClinGen Allele
Registry ID 識別變異；我們的 TSV 已有 CLINVAR_VARIATION_ID 欄位，直接對得上）。

輸入：ERepo 下載檔（TSV／CSV／JSON，可為 .gz）
      erepo.clinicalgenome.org → 搜尋頁或首頁的 "Download"（TSV summary），
      或 API：/evrepo/api/classifications
輸出：clingen_erepo_lookup.tsv.gz
      欄位：variation_id  class  criteria  panel  date

欄名採**模糊比對**（不寫死）：ERepo 匯出欄名各版本略有差異，故以關鍵字評分挑選對應欄，
並在 stderr 印出實際挑中的欄名供人工確認。

相依：只用 Python 標準庫。
"""

import argparse
import csv
import gzip
import io
import json
import sys


def _open_text(path):
    """支援 .gz 與純文字。"""
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _pick(headers, must_any, prefer_any=(), exclude_any=()):
    """
    從 headers 挑一個最像的欄名。
      must_any    ：欄名需包含其中任一關鍵字（缺一不可的條件）
      prefer_any  ：包含這些關鍵字者加分
      exclude_any ：包含這些關鍵字者直接排除
    回傳欄名或 None。
    """
    best, best_score = None, -1
    for h in headers:
        hl = (h or "").strip().lower()
        if not hl:
            continue
        if any(x in hl for x in exclude_any):
            continue
        if not any(x in hl for x in must_any):
            continue
        score = sum(2 for x in prefer_any if x in hl)
        score += 1 if len(hl) < 40 else 0      # 較短、直白的欄名優先
        if score > best_score:
            best, best_score = h, score
    return best


def detect_columns(headers):
    """回傳 {role: 實際欄名}，role = variation_id / class / criteria / panel / date。"""
    cols = {}
    # ClinVar Variation ID（排除 Allele Registry / VCI 之類的其他 ID）
    cols["variation_id"] = _pick(
        headers,
        must_any=("variation id", "variationid", "clinvar variation", "clinvar id",
                  "clinvar_variation_id", "variation"),
        prefer_any=("clinvar", "variation id"),
        exclude_any=("allele registry", "caid", "vci", "interpretation id"),
    )
    # 判讀結論
    cols["class"] = _pick(
        headers,
        must_any=("assertion", "classification", "pathogenicity", "significance"),
        prefer_any=("assertion", "classification"),
        exclude_any=("criteria", "code", "date", "count"),
    )
    # 套用的 ACMG criteria（要「Met」那一欄，排除 Not Met）
    cols["criteria"] = _pick(
        headers,
        must_any=("criteria", "evidence code", "evidence codes", "applied evidence"),
        prefer_any=("met",),
        exclude_any=("not met", "notmet"),
    )
    # VCEP / affiliation
    cols["panel"] = _pick(
        headers,
        must_any=("panel", "affiliation", "expert", "submitter", "vcep"),
        prefer_any=("panel", "vcep"),
        exclude_any=("id",),
    )
    # 判讀日期
    cols["date"] = _pick(
        headers,
        must_any=("date", "evaluated", "released"),
        prefer_any=("evaluated", "released"),
        exclude_any=(),
    )
    return cols


def _clean_varid(v):
    """把 '12345' / 'VCV000012345' / 'ClinVar:12345' 一律正規化成純數字字串。"""
    s = str(v or "").strip()
    if not s:
        return ""
    for pre in ("clinvar:", "vcv", "variation:"):
        if s.lower().startswith(pre):
            s = s[len(pre):]
    s = s.lstrip("0") or "0"
    return s if s.isdigit() else ""


def _clean(v):
    s = str(v or "").strip().replace("\t", " ").replace("\n", " ")
    return s if s else "."


def iter_records(path):
    """逐筆吐出 dict（自動判斷 JSON / TSV / CSV）。"""
    with _open_text(path) as f:
        head = f.read(4096)
        f.seek(0)
        stripped = head.lstrip()
        if stripped.startswith("{") or stripped.startswith("["):
            data = json.load(f)
            if isinstance(data, dict):
                for key in ("rows", "data", "classifications", "variantInterpretations", "items"):
                    if isinstance(data.get(key), list):
                        data = data[key]
                        break
                else:
                    data = [data]
            for rec in data:
                if isinstance(rec, dict):
                    yield _flatten(rec)
            return
        delim = "\t" if head.count("\t") >= head.count(",") else ","
        for rec in csv.DictReader(f, delimiter=delim):
            yield rec


def _flatten(d, prefix="", out=None):
    """JSON 巢狀結構壓平成 'a.b' 欄名，讓同一套欄名比對邏輯可共用。"""
    if out is None:
        out = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            _flatten(v, key + ".", out)
        elif isinstance(v, list):
            flat = [x for x in v if not isinstance(x, (dict, list))]
            if flat:
                out[key] = ", ".join(str(x) for x in flat)
        else:
            out[key] = v
    return out


def build(in_path, out_path):
    n_in = n_out = n_no_id = 0
    cols = None
    seen = {}

    for rec in iter_records(in_path):
        n_in += 1
        if cols is None:
            cols = detect_columns(list(rec.keys()))
            print("[build_clingen_erepo] 偵測到的欄位對應：", file=sys.stderr)
            for role in ("variation_id", "class", "criteria", "panel", "date"):
                print("    %-13s -> %s" % (role, cols.get(role) or "(找不到)"), file=sys.stderr)
            if not cols.get("variation_id") or not cols.get("class"):
                print("[ERROR] 找不到 ClinVar Variation ID 或 classification 欄位。"
                      "請把下載檔的 header 貼出來以便調整比對規則。\n"
                      "        現有欄名：%s" % list(rec.keys())[:40], file=sys.stderr)
                sys.exit(1)

        vid = _clean_varid(rec.get(cols["variation_id"]))
        if not vid:
            n_no_id += 1
            continue

        row = (
            _clean(rec.get(cols["class"])),
            _clean(rec.get(cols["criteria"])) if cols.get("criteria") else ".",
            _clean(rec.get(cols["panel"])) if cols.get("panel") else ".",
            _clean(rec.get(cols["date"])) if cols.get("date") else ".",
        )
        # 同一變異可能有多個 VCEP 判讀：保留第一筆，另記錄有幾筆（供人工複核）
        if vid in seen:
            prev = seen[vid]
            seen[vid] = (prev[0], prev[1], prev[2], prev[3], prev[4] + 1)
        else:
            seen[vid] = (row[0], row[1], row[2], row[3], 1)
            n_out += 1

    with gzip.open(out_path, "wt") as w:
        w.write("variation_id\tclass\tcriteria\tpanel\tdate\tn_records\n")
        for vid in sorted(seen, key=lambda x: int(x)):
            c, crit, panel, date, n = seen[vid]
            w.write("%s\t%s\t%s\t%s\t%s\t%d\n" % (vid, c, crit, panel, date, n))

    print("[build_clingen_erepo] 讀入 %d 筆，輸出 %d 個變異（無 ClinVar ID 略過 %d 筆）"
          % (n_in, n_out, n_no_id), file=sys.stderr)
    multi = sum(1 for v in seen.values() if v[4] > 1)
    if multi:
        print("[build_clingen_erepo] 其中 %d 個變異有多個 VCEP 判讀（只保留第一筆，"
              "n_records 欄記錄筆數）" % multi, file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(
        description="Convert a ClinGen Evidence Repository download into a Variation-ID keyed lookup.")
    ap.add_argument("--input", required=True, help="ERepo download (TSV/CSV/JSON, may be .gz)")
    ap.add_argument("--output", required=True, help="clingen_erepo_lookup.tsv.gz")
    a = ap.parse_args()
    build(a.input, a.output)


if __name__ == "__main__":
    main()
