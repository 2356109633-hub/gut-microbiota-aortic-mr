# ============================================================================
# 肠道菌群与主动脉疾病双向 MR —— FinnGen R13(DF13) 发现 + AAAgen 验证（双队列）
#   v2 修正版   2026-10-25
#
# 相对 v1 (mr_R13_AAAgen_dual.R) 的四处修正 —— 详见《验证报告_R13_AAAgen_独立复算.md》
#   ① 暴露端剔除伪暴露 Reverse_AD（原文件混入 95 行 / 32 个 SNP 的 AD 结局数据）
#   ② AAAgen 放弃手写位置匹配：用暴露 hg19 坐标给 AAAgen 变异赋伪 rsID，交给
#      TwoSampleMR::harmonise_data() 按等位基因方向自动翻号。
#      v1 的 harmonise_by_pos() 在等位基因互换时不翻号 → 约 47% 的 IV 以反号进入 IVW，
#      导致 v1 的 AAAgen 结果整体失效。
#   ③ 敏感性分析改为对全部 taxa 执行（v1 挂在"双队列交集"上，交集为空即整段跳过）
#   ④ FDR 按 140 个 taxa 计算；dual 表无论是否为空都写出（避免残留旧文件）
#
# 环境：R 4.3.1 + TwoSampleMR 0.4.26（本机实际版本）
# 用法：RStudio 打开本文件 → Session → Set Working Directory → To Source File Location → Source
#      预计 10~20 分钟（主要耗时在读 4 个 FinnGen gz 与 AAAgen 全表）
# ============================================================================

## ---------------------- 0. 配置 ----------------------
CLEAR_ENV <- TRUE
if (CLEAR_ENV) { rm(list = ls(all.names = TRUE)); gc(); cat("[INFO] 已清空 Global Environment\n") }

outdir <- "results_R13_AAAgen_v2"
dir.create(outdir, showWarnings = FALSE)
LOG_FILE <- file.path(outdir, "run_log_R13_AAAgen_v2.txt")
if (file.exists(LOG_FILE)) file.remove(LOG_FILE)
sink(LOG_FILE, split = TRUE)
# 注意：顶层 on.exit() 不会执行，日志必须在脚本末尾显式关闭（见文件最后）。
# 若在 RStudio 里中途按 Esc 中断，请手动运行：while (sink.number() > 0) sink()

cat("[INFO] 运行时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("[INFO] 输出目录:", outdir, "\n")

fg_dir <- "finngen-public-data-r13"
aaagen_dir <- c("AAAGen", "AAgen", "AAAgen")[file.exists(c("AAAGen", "AAgen", "AAAgen"))][1]
if (is.na(aaagen_dir)) stop("未找到 AAAgen 数据文件夹（AAAGen / AAgen / AAAgen）")

# 病例/对照数（FinnGen R13 官方 manifest 核实；AAAgen 取 Roychowdhury 2023）
case_info <- list(
  AAA_R13    = c(4634,   460490),
  TAA_R13    = c(5704,   460490),
  AA_R13     = c(10694,  460490),
  AD_R13     = c(1219,   460490),
  AAA_AAAgen = c(39221, 1086107)   # 39,221 cases total (37,214 EUR + 2,007 AFR), per Roychowdhury 2023
)
N_EXP <- 18340   # MiBioGen

## ---------------------- 1. 包 ----------------------
suppressMessages({ library(TwoSampleMR); library(data.table); library(dplyr) })
set.seed(2026)
cat("[INFO] TwoSampleMR", as.character(packageVersion("TwoSampleMR")),
    "| data.table", as.character(packageVersion("data.table")), "\n")
cat("[INFO]", R.version.string, "\n")

## ---------------------- 2. 暴露 ----------------------
NON_TAXA <- c("Reverse_AD")   # ① 伪暴露

load_exposure <- function(f = "harmonised_dat.csv") {
  if (!file.exists(f)) stop("未找到 ", f)
  d <- fread(f)
  setnames(d, names(d), tolower(names(d)))
  setnames(d, names(d), sub("\\.exposure$", "", names(d)))
  if ("snp" %in% names(d) && !"SNP" %in% names(d)) setnames(d, "snp", "SNP")
  if (!"eaf" %in% names(d)) { d[, eaf := NA_real_]; cat("[INFO] 暴露表无 eaf，已置 NA（不影响主分析）\n") }
  d <- d[!is.na(SNP) & !is.na(beta) & !is.na(se)]
  d <- d[order(exposure, SNP, pval)]
  d <- d[!duplicated(paste(exposure, SNP))]
  n0 <- uniqueN(d$exposure)
  d <- d[!exposure %in% NON_TAXA]
  cat("[OK] 暴露: 剔除伪暴露", paste(NON_TAXA, collapse = ","), "->",
      n0, "→", uniqueN(d$exposure), "个 taxa\n")

  # hg19 坐标（AAAgen 为 hg19，必须用 hg19；hg38 会全部失配）
  if (!all(c("chrom", "pos") %in% names(d)) || all(is.na(d$pos))) {
    if (all(c("chrom", "pos") %in% names(d))) d[, c("chrom", "pos") := NULL]
    cf <- "snps_hg19.csv"
    if (!file.exists(cf) || file.info(cf)$size == 0) stop("缺少 snps_hg19.csv（hg19 坐标）")
    b <- fread(cf)
    snp_col <- which(grepl("snp|rsid|^rs$|variant", names(b), ignore.case = TRUE))[1]
    chr_col <- which(grepl("^chr|chrom", names(b), ignore.case = TRUE))[1]
    pos_col <- which(grepl("pos", names(b), ignore.case = TRUE) &
                       !grepl("snp|rs|variant", names(b), ignore.case = TRUE))[1]
    b2 <- data.table(SNP   = as.character(b[[snp_col]]),
                     chrom = sub("^chr", "", as.character(b[[chr_col]])),
                     pos   = suppressWarnings(as.numeric(b[[pos_col]])))
    b2 <- b2[!is.na(SNP) & !is.na(pos)][!duplicated(SNP)]
    d <- merge(d, b2, by = "SNP", all.x = TRUE)
    cat("[INFO] hg19 坐标补充:", sum(!is.na(d$pos)), "/", nrow(d), "个 taxon-SNP 对\n")
  }
  d
}

exp_raw <- load_exposure()
exp_dat <- format_data(
  as.data.frame(exp_raw), type = "exposure", phenotype_col = "exposure", snp_col = "SNP",
  beta_col = "beta", se_col = "se", eaf_col = "eaf",
  effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "pval")
exp_dat$samplesize.exposure <- N_EXP
exp_dat$units.exposure <- "SD"

exp_dt   <- as.data.table(exp_dat)
iv_count <- exp_dt[, .N, by = exposure]
cat("[OK] 暴露 IV:", nrow(exp_dt), "个 taxon-SNP 对 /", uniqueN(exp_dt$SNP), "个唯一 SNP /",
    uniqueN(exp_dt$exposure), "个 taxa\n")
cat("[OK] 每 taxa IV 数: 中位数", median(iv_count$N), "，范围", min(iv_count$N), "-", max(iv_count$N), "\n")

# 位置 -> rsID（供 AAAgen 赋伪 rsID）
pos_tbl <- unique(exp_raw[!is.na(pos), .(SNP, chrom = as.character(chrom), pos)])
cat("[OK] 有 hg19 坐标的 SNP:", uniqueN(pos_tbl$SNP), "/", uniqueN(exp_raw$SNP), "\n")

## ---------------------- 3. FinnGen R13 结局 ----------------------
scan_fg <- function(dir) {
  files <- list.files(dir, pattern = "finngen_R13_I9_.*\\.gz$", full.names = TRUE)
  if (length(files) == 0) stop("未在 ", dir, " 找到 finngen_R13_I9_*.gz")
  map <- c(ABAORTANEUR = "AAA_R13", THAORTANEUR = "TAA_R13",
           AORTANEUR = "AA_R13", AORTDIS = "AD_R13")
  out <- list()
  for (pat in names(map)) {
    f <- files[grepl(paste0("I9_", pat, "\\.gz$"), files)]
    if (length(f) > 0) out[[map[pat]]] <- f[1]
  }
  out
}
fg_files <- scan_fg(fg_dir)
cat("[OK] FinnGen 文件:", paste(names(fg_files), collapse = ", "), "\n")

fmt_fg <- function(f, nm) {
  sel <- c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt")
  d <- fread(f, select = sel)
  if ("#chrom" %in% names(d)) setnames(d, "#chrom", "chrom")
  d <- d[rsids %in% unique(exp_dt$SNP)]
  cat("[INFO]", nm, "预过滤后:", nrow(d), "行\n")
  d[, outcome := nm][, SNP := rsids][, effect_allele := alt][, other_allele := ref]
  d[, eaf := af_alt][, se := sebeta]
  d <- d[!is.na(SNP) & !is.na(beta) & !is.na(se) & !is.na(eaf)]
  o <- format_data(as.data.frame(d), type = "outcome", phenotype_col = "outcome", snp_col = "SNP",
                   beta_col = "beta", se_col = "se", eaf_col = "eaf",
                   effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                   pval_col = "pval")
  o$ncase.outcome <- case_info[[nm]][1]
  o$ncontrol.outcome <- case_info[[nm]][2]
  o$samplesize.outcome <- sum(case_info[[nm]])
  o$units.outcome <- "log odds"
  o
}
fg_outcomes <- rbindlist(lapply(names(fg_files), function(nm) fmt_fg(fg_files[[nm]], nm)), fill = TRUE)
cat("[OK] FinnGen 结局合并:", nrow(fg_outcomes), "行 /",
    uniqueN(fg_outcomes$SNP), "SNP /", uniqueN(fg_outcomes$outcome), "个结局\n")

## ---------------------- 4. AAAgen（② 伪 rsID 法） ----------------------
aaagen_file <- list.files(aaagen_dir, pattern = "\\.(txt|gz|tsv)$", full.names = TRUE)[1]
cat("[INFO] AAAgen 文件:", aaagen_file, "\n")
probe <- fread(aaagen_file, nrows = 5)
need <- c("chr", "pos", "ref", "alt", "AF", "Effectsize", "Effectsize_SD", "pvalue")
if (!all(need %in% names(probe))) {
  stop("AAAgen 列名不符，实际为: ", paste(names(probe), collapse = ", "))
}
ag <- fread(aaagen_file, select = need)
cat("[OK] AAAgen 全表读入:", nrow(ag), "行\n")
gc()

ag[, chrom := sub("^chr", "", as.character(chr))]
ag[, pos := suppressWarnings(as.numeric(pos))]
ag[, effect_allele := toupper(alt)][, other_allele := toupper(ref)]

pm <- unique(pos_tbl[, .(chrom, pos, SNP)])
pos_set <- unique(pm$pos)
cat("[DIAG] 坐标表:", nrow(pm), "个 (chrom,pos)->SNP；唯一位置:", length(pos_set), "\n")

n_before <- nrow(ag)
ag <- ag[pos %in% pos_set]                      # 数值比较，避免字符串键的格式陷阱
cat("[INFO] 按 hg19 位置预过滤:", n_before, "→", nrow(ag), "行\n")
if (nrow(ag) == 0) stop("AAAgen 与暴露位置零匹配：坐标版本可能不是 hg19")

ag <- merge(ag, pm, by = c("chrom", "pos"))     # 两列连接，取回伪 rsID
cat("[DIAG] 两列连接后:", nrow(ag), "行；唯一 SNP:", uniqueN(ag$SNP),
    "；唯一位置:", uniqueN(paste(ag$chrom, ag$pos)), "\n")
if (nrow(ag) == 0) stop("位置命中后连接失败")

# 同一位置多行（多等位/重复记录）→ 优先保留与暴露等位基因对得上的那一行
ea <- unique(exp_raw[, .(SNP, ea_exp = toupper(effect_allele), oa_exp = toupper(other_allele))])
ag <- merge(ag, ea, by = "SNP", all.x = TRUE)
ag[, match_rank := fifelse(effect_allele == ea_exp & other_allele == oa_exp, 1L,
                    fifelse(effect_allele == oa_exp & other_allele == ea_exp, 2L, 3L))]
setorder(ag, SNP, match_rank, pvalue)
n_multi <- ag[, .N, by = SNP][N > 1, .N]
ag <- ag[!duplicated(SNP)]
cat("[INFO] 一位点多记录:", n_multi, "个位置；保留结果 —— 同向", ag[match_rank == 1, .N],
    "/ 互换", ag[match_rank == 2, .N], "/ 不匹配", ag[match_rank == 3, .N], "\n")

ag[, eaf := AF][, beta := Effectsize][, se := Effectsize_SD][, pval := pvalue][, outcome := "AAA_AAAgen"]
aaagen_dat <- format_data(as.data.frame(ag), type = "outcome", phenotype_col = "outcome", snp_col = "SNP",
                          beta_col = "beta", se_col = "se", eaf_col = "eaf",
                          effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                          pval_col = "pval")
aaagen_dat$ncase.outcome <- case_info[["AAA_AAAgen"]][1]
aaagen_dat$ncontrol.outcome <- case_info[["AAA_AAAgen"]][2]
aaagen_dat$samplesize.outcome <- sum(case_info[["AAA_AAAgen"]])
aaagen_dat$units.outcome <- "log odds"
cat("[OK] AAAgen 可用变异:", nrow(aaagen_dat), "\n")
rm(ag); gc()

## ---------------------- 5. 工具函数 ----------------------
# 注意：TwoSampleMR 的 harmonise_data / mr 会打印大量进度信息，全部静音
MT <- c("mr_ivw", "mr_weighted_median", "mr_egger_regression")

harm <- function(ed, od) {
  invisible(capture.output(
    dat <- tryCatch(harmonise_data(as.data.frame(ed), as.data.frame(od)), error = function(e) NULL)))
  if (is.null(dat) || !"mr_keep" %in% names(dat)) return(NULL)
  d <- as.data.table(dat)[mr_keep == TRUE]
  if (nrow(d) < 3) return(NULL)
  d
}

do_mr <- function(d) {
  invisible(capture.output(
    r <- tryCatch(as.data.table(mr(as.data.frame(d), method_list = MT)), error = function(e) NULL)))
  if (is.null(r)) return(NULL)
  r[, n_used := nrow(d)]
  r
}

## ---------------------- 6. 正向 MR：R13 四结局 ----------------------
cat("\n[RUN] 正向 R13:", uniqueN(exp_dt$exposure), "taxa x", length(fg_files), "结局\n")
res13 <- list()
for (oc in names(fg_files)) {
  od <- fg_outcomes[outcome == oc]
  for (tx in unique(exp_dt$exposure)) {
    d <- harm(exp_dt[exposure == tx], od)
    if (is.null(d)) next
    r <- do_mr(d)
    if (is.null(r)) next
    r[, `:=`(outcome = oc, exposure = tx)]
    res13[[paste(tx, oc)]] <- r
  }
  cat("  [", oc, "] 累计", length(res13), "条\n")
}
fw13 <- rbindlist(res13, fill = TRUE)
cat("[OK] R13 正向:", nrow(fw13), "行 /", uniqueN(fw13$exposure), "taxa\n")
write.csv(fw13, file.path(outdir, "forward_R13_all.csv"), row.names = FALSE)

## ---------------------- 7. 正向 MR：AAAgen ----------------------
cat("\n[RUN] 正向 AAAgen\n")
resA <- list()
for (tx in unique(exp_dt$exposure)) {
  d <- harm(exp_dt[exposure == tx], aaagen_dat)
  if (is.null(d)) next
  r <- do_mr(d)
  if (is.null(r)) next
  r[, `:=`(outcome = "AAA_AAAgen", exposure = tx)]
  resA[[tx]] <- r
}
fwA <- rbindlist(resA, fill = TRUE)
cat("[OK] AAAgen 正向:", nrow(fwA), "行 /", uniqueN(fwA$exposure), "taxa\n")
write.csv(fwA, file.path(outdir, "forward_AAAgen_all.csv"), row.names = FALSE)

## ---------------------- 8. FDR / Bonferroni 表 ----------------------
add_fdr <- function(dt, by_outcome) {
  ivw <- dt[method == "Inverse variance weighted"]
  if (by_outcome) {
    ivw[, `:=`(q_fdr = p.adjust(pval, "fdr"), q_bonf = p.adjust(pval, "bonferroni")), by = outcome]
  } else {
    ivw[, `:=`(q_fdr = p.adjust(pval, "fdr"), q_bonf = p.adjust(pval, "bonferroni"))]
  }
  ivw[, `:=`(OR = exp(b), ci_low = exp(b - 1.96 * se), ci_high = exp(b + 1.96 * se))]
  ivw
}
fdr13 <- add_fdr(fw13, TRUE)
fdrA  <- add_fdr(fwA,  FALSE)
write.csv(fdr13, file.path(outdir, "forward_R13_ivw_fdr.csv"), row.names = FALSE)
write.csv(fdrA,  file.path(outdir, "forward_AAAgen_ivw_fdr.csv"), row.names = FALSE)

cat("\n[摘要] IVW 名义显著 taxa 数 ——\n")
print(fdr13[pval < 0.05, .N, by = outcome])
cat("[摘要] AAAgen:", fdrA[pval < 0.05, .N], "个\n")
cat("[摘要] 最小 q_fdr —— R13:", signif(min(fdr13$q_fdr), 4),
    "| AAAgen:", signif(min(fdrA$q_fdr), 4), "\n")

## ---------------------- 9. 双队列交叉核对 ----------------------
cat("\n[RUN] 双队列交叉核对\n")
a <- fdr13[outcome == "AAA_R13", .(exposure, n_R13 = nsnp, b_R13 = b, se_R13 = se,
                                   p_R13 = pval, q_R13 = q_fdr)]
b <- fdrA[, .(exposure, n_AG = nsnp, b_AG = b, se_AG = se, p_AG = pval, q_AG = q_fdr)]
xc <- merge(a, b, by = "exposure", all = TRUE)
xc[, same_dir := sign(b_R13) == sign(b_AG)]
xc[, `:=`(OR_R13 = exp(b_R13), OR_AG = exp(b_AG),
          ci_R13_low = exp(b_R13 - 1.96 * se_R13), ci_R13_high = exp(b_R13 + 1.96 * se_R13),
          ci_AG_low  = exp(b_AG  - 1.96 * se_AG),  ci_AG_high  = exp(b_AG  + 1.96 * se_AG))]
xc[, robust := (q_R13 < 0.05) & (q_AG < 0.05) & same_dir]
setorder(xc, p_R13)
write.csv(xc, file.path(outdir, "cross_cohort_lookup.csv"), row.names = FALSE)

# ④ dual 表无条件写出（空表也要有表头，避免误读旧文件）
dual <- xc[!is.na(p_R13) & !is.na(p_AG) & p_R13 < 0.05 & p_AG < 0.05 & same_dir == TRUE]
write.csv(dual, file.path(outdir, "dual_validation_AAA.csv"), row.names = FALSE)
cat("[OK] 共同 taxa:", nrow(xc),
    "| 方向一致:", xc[same_dir == TRUE, .N],
    "(", round(100 * xc[same_dir == TRUE, .N] / xc[!is.na(same_dir), .N], 1), "%)",
    "| 双侧名义显著:", nrow(dual),
    "| FDR 稳健:", xc[robust == TRUE, .N], "\n")

## ---------------------- 10. 敏感性分析（全部 taxa） ----------------------
cat("\n[RUN] 敏感性分析：全部 taxa x 5 结局\n")
sens_one <- function(ed, od, tx, oc) {
  d <- harm(ed, od)
  if (is.null(d)) return(NULL)
  invisible(capture.output(
    het <- tryCatch(as.data.table(mr_heterogeneity(as.data.frame(d), method_list = "mr_ivw")),
                    error = function(e) NULL)))
  invisible(capture.output(
    ple <- tryCatch(as.data.table(mr_pleiotropy_test(as.data.frame(d))), error = function(e) NULL)))
  tstat <- d$beta.exposure / d$se.exposure
  r2e <- sum(tstat^2 / (tstat^2 + N_EXP - 2))          # 连续暴露由 t 反推 r
  prev <- d$ncase.outcome[1] / (d$ncase.outcome[1] + d$ncontrol.outcome[1])
  r2o <- tryCatch({
    r <- get_r_from_lor(d$beta.outcome, d$eaf.outcome, d$ncase.outcome, d$ncontrol.outcome, prev)
    sum(r^2) }, error = function(e) NA_real_)
  data.table(exposure = tx, outcome = oc, nsnp = nrow(d),
             Q = if (!is.null(het)) het$Q[1] else NA_real_,
             Q_df = if (!is.null(het)) het$Q_df[1] else NA_real_,
             Q_pval = if (!is.null(het)) het$Q_pval[1] else NA_real_,
             Egger_intercept = if (!is.null(ple)) ple$egger_intercept[1] else NA_real_,
             Egger_se = if (!is.null(ple)) ple$se[1] else NA_real_,
             Egger_pval = if (!is.null(ple)) ple$pval[1] else NA_real_,
             R2_exposure = r2e, R2_outcome = r2o,
             steiger_dir = ifelse(!is.na(r2o), r2e > r2o, NA))
}
sens <- list()
for (oc in names(fg_files)) {
  od <- fg_outcomes[outcome == oc]
  for (tx in unique(exp_dt$exposure)) {
    s <- sens_one(exp_dt[exposure == tx], od, tx, oc)
    if (!is.null(s)) sens[[paste(tx, oc)]] <- s
  }
  cat("  [", oc, "] 完成\n")
}
for (tx in unique(exp_dt$exposure)) {
  s <- sens_one(exp_dt[exposure == tx], aaagen_dat, tx, "AAA_AAAgen")
  if (!is.null(s)) sens[[paste(tx, "AG")]] <- s
}
sens <- rbindlist(sens, fill = TRUE)
write.csv(sens, file.path(outdir, "sensitivity_all.csv"), row.names = FALSE)

cand <- unique(c(fdr13[pval < 0.05, exposure], fdrA[pval < 0.05, exposure]))
write.csv(sens[exposure %in% cand], file.path(outdir, "sensitivity_candidates.csv"), row.names = FALSE)
cat("[OK] 敏感性:", nrow(sens), "行；候选 taxa", length(cand), "个\n")

## ---------------------- 10b. 留一法 与 MR-PRESSO（候选 taxa） ----------------------
cat("\n[RUN] 留一法 (leave-one-out) 与 MR-PRESSO：候选 taxa\n")
# MRPRESSO 为可选依赖：未安装则跳过并记录，不中断主流程。
# 安装方式（需联网，在本机 R 控制台执行一次即可）：install.packages("MRPRESSO")
has_presso <- requireNamespace("MRPRESSO", quietly = TRUE)
cat("[INFO] MRPRESSO 可用:", has_presso,
    if (!has_presso) "（未安装 → 跳过 MR-PRESSO，其余检验不受影响）" else "", "\n")

loo_list <- list(); presso_list <- list()
sens_outcomes <- c(names(fg_files), "AAA_AAAgen")
for (oc in sens_outcomes) {
  od <- if (oc == "AAA_AAAgen") aaagen_dat else fg_outcomes[outcome == oc]
  for (tx in cand) {
    d <- harm(exp_dt[exposure == tx], od)
    if (is.null(d)) next

    # --- 留一法 ---
    invisible(capture.output(
      lo <- tryCatch(as.data.table(mr_leaveoneout(as.data.frame(d))), error = function(e) NULL)))
    if (!is.null(lo) && nrow(lo) > 0) {
      lo[, `:=`(exposure = tx, outcome = oc)]
      loo_list[[paste(tx, oc)]] <- lo
    }

    # --- MR-PRESSO（≥4 个 SNP 才有意义） ---
    if (has_presso && nrow(d) >= 4) {
      set.seed(2026)
      # 注意：不要用 capture.output() 包裹 mr_presso —— 其内部会接管输出，
      # 嵌套 sink 会导致报错并静默丢结果。这里只用 suppressWarnings/suppressMessages。
      pr <- suppressWarnings(suppressMessages(tryCatch(MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome", SdExposure = "se.exposure",
        data = as.data.frame(d), OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
        NbDistribution = 1000, SignifThreshold = 0.05), error = function(e) NULL)))
      if (!is.null(pr)) {
        # MRPRESSO 封装：[[1]] Main MR results (data.frame) / [[2]] MR-PRESSO results (list)
        # 未检出离群点 时 list 中只含 Global Test（无 Outlier/Distortion Test）
        pres <- tryCatch(pr[["MR-PRESSO results"]], error = function(e) NULL)
        g    <- tryCatch(pres[["Global Test"]]$Pvalue, error = function(e) NA_real_)
        ot   <- tryCatch(pres[["Outlier Test"]], error = function(e) NULL)
        n_out <- if (!is.null(ot)) sum(ot$Pvalue < 0.05, na.rm = TRUE) else 0L
        out_snps <- if (!is.null(ot)) paste(rownames(ot)[ot$Pvalue < 0.05], collapse = ";") else ""
        dt <- tryCatch(pres[["Distortion Test"]]$Pvalue, error = function(e) NA_real_)
        mr_tab <- tryCatch(as.data.table(pr[["Main MR results"]]), error = function(e) NULL)
        b_raw <- if (!is.null(mr_tab)) mr_tab[`MR Analysis` == "Raw", `Causal Estimate`][1] else NA_real_
        p_raw <- if (!is.null(mr_tab)) mr_tab[`MR Analysis` == "Raw", `P-value`][1] else NA_real_
        b_cor <- if (!is.null(mr_tab)) mr_tab[`MR Analysis` == "Outlier-corrected", `Causal Estimate`][1] else NA_real_
        p_cor <- if (!is.null(mr_tab)) mr_tab[`MR Analysis` == "Outlier-corrected", `P-value`][1] else NA_real_
        presso_list[[paste(tx, oc)]] <- data.table(
          exposure = tx, outcome = oc, nsnp = nrow(d),
          presso_global_p = as.numeric(g),
          n_outliers = n_out,
          outlier_snps = out_snps,
          distortion_p = if (is.null(dt)) NA_real_ else as.numeric(dt)[1],
          b_raw = as.numeric(b_raw), p_raw = as.numeric(p_raw),
          b_outlier_corrected = as.numeric(b_cor), p_outlier_corrected = as.numeric(p_cor))
      }
    }
  }
  cat("  [", oc, "] LOO/PRESSO 完成\n")
}

if (length(loo_list) > 0) {
  loo <- rbindlist(loo_list, fill = TRUE)
  write.csv(loo, file.path(outdir, "loo_candidates.csv"), row.names = FALSE)
  # 与全模型比较：是否有单个 SNP 主导结果
  ivw_full <- rbind(
    fdr13[exposure %in% cand, .(exposure, outcome, b_all = b, p_all = pval)],
    fdrA[exposure %in% cand, .(exposure, outcome = "AAA_AAAgen", b_all = b, p_all = pval)]
  )
  # 注意：mr_leaveoneout() 的 P 值列名是 `p`（不是 pval）
  loo_sum <- loo[, .(n_loo = .N, b_loo_min = min(b), b_loo_max = max(b),
                     p_loo_max = max(p)), by = .(exposure, outcome)]
  loo_sum <- merge(loo_sum, ivw_full, by = c("exposure", "outcome"), all.x = TRUE)
  loo_sum[, max_abs_delta := pmax(abs(b_loo_min - b_all), abs(b_loo_max - b_all))]
  loo_sum[, driven_by_single_snp := (p_all < 0.05) & (p_loo_max > 0.05)]
  setorder(loo_sum, -max_abs_delta)
  write.csv(loo_sum, file.path(outdir, "loo_summary.csv"), row.names = FALSE)
  cat("[OK] 留一法:", nrow(loo), "行；汇总:", nrow(loo_sum), "行\n")
  cat("[INFO] 存在单 SNP 即可推翻显著性的候选:", loo_sum[driven_by_single_snp == TRUE, .N], "个\n")
} else {
  cat("[!] 留一法无结果\n")
}

if (length(presso_list) > 0) {
  presso <- rbindlist(presso_list, fill = TRUE)
  write.csv(presso, file.path(outdir, "presso_candidates.csv"), row.names = FALSE)
  cat("[OK] MR-PRESSO:", nrow(presso), "行；全局检验 P<0.05:",
      presso[!is.na(presso_global_p) & presso_global_p < 0.05, .N], "个\n")
} else {
  cat("[!] MR-PRESSO 未运行（未安装 MRPRESSO 或所有 taxa SNP<4），未生成 presso_candidates.csv\n")
}

## ---------------------- 11. 反向 MR ----------------------
# 限制说明：疾病 IV 只能从 MiBioGen IV 池中挑（本脚本只有该池的 SNP 数据），
# 且 taxon 作结局时仅有其自身 IV 可用 → 覆盖不完整。规范反向 MR 需 MiBioGen 全量
# summary statistics。论文 Limitations 必须如实说明这一点。
cat("\n[RUN] 反向 MR（受限）\n")
extract_rev <- function(od, nm, pthr = 5e-6) {
  iv <- od[pval.outcome < pthr & !is.na(beta.outcome) & !is.na(se.outcome)]
  if (nrow(iv) == 0) return(NULL)
  iv[, .(SNP, effect_allele = effect_allele.outcome, other_allele = other_allele.outcome,
         eaf = eaf.outcome, beta = beta.outcome, se = se.outcome, pval = pval.outcome,
         exposure = nm)]
}
rev_raw <- rbindlist(lapply(names(fg_files), function(nm) extract_rev(fg_outcomes[outcome == nm], nm)),
                     fill = TRUE)
if (is.null(rev_raw)) rev_raw <- data.table()
if (nrow(rev_raw) > 0) rev_raw <- rev_raw[!duplicated(paste(exposure, SNP))]
cat("[INFO] 各结局可用疾病 IV 数:\n")
if (nrow(rev_raw) > 0) print(rev_raw[, .N, by = exposure])

if (nrow(rev_raw) > 0) {
  rev_exp <- format_data(as.data.frame(rev_raw), type = "exposure", phenotype_col = "exposure",
                         snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf",
                         effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                         pval_col = "pval")
  rev_out_raw <- copy(exp_raw)[, outcome := paste0("taxa_", exposure)]   # exp_raw 已剔除 Reverse_AD
  rev_out <- format_data(as.data.frame(rev_out_raw), type = "outcome", phenotype_col = "outcome",
                         snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf",
                         effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                         pval_col = "pval")
  rl <- list()
  for (ex in unique(rev_exp$exposure)) {
    ed <- as.data.table(rev_exp)[exposure == ex]
    for (oc in unique(rev_out$outcome)) {
      od <- as.data.table(rev_out)[outcome == oc]
      d <- harm(ed, od)
      if (is.null(d)) next
      r <- do_mr(d)
      if (is.null(r)) next
      r <- r[method == "Inverse variance weighted"]
      if (nrow(r) == 0) next
      r[, `:=`(exposure = ex, outcome = oc, n_used = d[, .N])]
      rl[[paste(ex, oc)]] <- r
    }
  }
  rev <- rbindlist(rl, fill = TRUE)
  write.csv(rev, file.path(outdir, "reverse_all.csv"), row.names = FALSE)
  cat("[OK] 反向 MR:", nrow(rev), "行\n")
  if (nrow(rev) > 0) {
    cat("[INFO] 覆盖结局:", paste(unique(rev$exposure), collapse = ", "), "\n")
    cat("[INFO] 名义 P<0.05:", rev[pval < 0.05, .N], "| 最小 P:", signif(min(rev$pval), 4), "\n")
    cat("[INFO] 未覆盖的结局:", paste(setdiff(names(fg_files), unique(rev$exposure)), collapse = ", "), "\n")
  }
} else {
  cat("[!] 无疾病 IV，跳过反向 MR\n")
}

## ---------------------- 12. 论文数字汇总 ----------------------
sumf <- file.path(outdir, "论文数字汇总_v2.txt")
if (file.exists(sumf)) file.remove(sumf)
sm <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  write(txt, file = sumf, append = TRUE)
}
sm("=== 双队列 MR 最终结果（v2 修正版） ===")
sm("日期: ", format(Sys.time(), "%Y-%m-%d %H:%M"))
sm("TwoSampleMR ", as.character(packageVersion("TwoSampleMR")), " + ", R.version.string)
sm("")
sm("--- 暴露 ---")
sm("MiBioGen: ", uniqueN(exp_dt$exposure), " taxa / ", nrow(exp_dt), " 个 taxon-SNP 对 / ",
   uniqueN(exp_dt$SNP), " 个唯一 SNP（每 taxa 中位 ", median(iv_count$N), " 个 IV）")
sm("")
sm("--- 结局 ---")
for (k in names(case_info)) {
  sm(k, ": cases=", format(case_info[[k]][1], big.mark = ","),
     ", controls=", format(case_info[[k]][2], big.mark = ","))
}
sm("")
sm("--- 各结局 IVW 名义显著 (P<0.05) ---")
for (oc in names(fg_files)) {
  s <- fdr13[outcome == oc & pval < 0.05][order(pval)]
  sm("[", oc, "] n=", nrow(s))
  if (nrow(s)) for (i in seq_len(nrow(s))) {
    sm("  ", s$exposure[i], ": b=", round(s$b[i], 4), ", se=", round(s$se[i], 4),
       ", p=", signif(s$pval[i], 4), ", OR=", round(s$OR[i], 3),
       " (", round(s$ci_low[i], 3), "-", round(s$ci_high[i], 3), "), q_fdr=", signif(s$q_fdr[i], 3))
  }
}
s <- fdrA[pval < 0.05][order(pval)]
sm("[AAAgen] n=", nrow(s))
if (nrow(s)) for (i in seq_len(nrow(s))) {
  sm("  ", s$exposure[i], ": b=", round(s$b[i], 4), ", se=", round(s$se[i], 4),
     ", p=", signif(s$pval[i], 4), ", OR=", round(s$OR[i], 3),
     " (", round(s$ci_low[i], 3), "-", round(s$ci_high[i], 3), "), q_fdr=", signif(s$q_fdr[i], 3))
}
sm("")
sm("--- 双队列交叉 ---")
sm("共同 taxa: ", nrow(xc), "；方向一致: ", xc[same_dir == TRUE, .N], "/", nrow(xc))
sm("双侧名义显著: ", nrow(dual), "；FDR 稳健: ", xc[robust == TRUE, .N])
sm("")
sm("--- 敏感性（候选 taxa） ---")
sc <- sens[exposure %in% cand][order(exposure, outcome)]
if (nrow(sc)) for (i in seq_len(nrow(sc))) {
  sm("  ", sc$exposure[i], " [", sc$outcome[i], "] Q=", round(sc$Q[i], 2),
     " (p=", signif(sc$Q_pval[i], 3), "), Egger截距=", round(sc$Egger_intercept[i], 4),
     " (p=", signif(sc$Egger_pval[i], 3), "), Steiger=",
     ifelse(isTRUE(sc$steiger_dir[i]), "方向正确", "方向异常"))
}
if (file.exists(file.path(outdir, "loo_summary.csv"))) {
  ls <- fread(file.path(outdir, "loo_summary.csv"))
  sm("")
  sm("--- 留一法（主结局：R13 AAA 与 AAAgen AAA） ---")
  ls2 <- ls[outcome %in% c("AAA_R13", "AAA_AAAgen")][order(outcome, p_all)]
  for (i in seq_len(nrow(ls2))) {
    sm("  ", ls2$exposure[i], " [", ls2$outcome[i], "] b=", round(ls2$b_all[i], 3),
       ", p=", signif(ls2$p_all[i], 3), ", 单 SNP 最大影响=", round(ls2$max_abs_delta[i], 3),
       ", LOO 最大 P=", signif(ls2$p_loo_max[i], 3),
       if (isTRUE(ls2$driven_by_single_snp[i])) " <- 单个 SNP 即可推翻" else " (稳定)")
  }
  sm("")
  sm("--- MR-PRESSO ---")
  sm(if (has_presso) "已运行（见 presso_candidates.csv）" else "未运行：本机未安装 MRPRESSO 包（需联网 install.packages(\"MRPRESSO\")）")
}

cat("\n[OK] 汇总已写入:", sumf, "\n")

cat("\n========== 完成 ==========\n")
cat("输出目录:", normalizePath(outdir), "\n")

# 显式关闭日志 sink（顶层 on.exit 不生效，必须在此关闭，否则 RStudio 后续输出被吞）
while (sink.number() > 0) sink()
cat("[OK] 日志已关闭:", normalizePath(LOG_FILE), "\n")
