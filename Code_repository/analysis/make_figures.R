# ============================================================================
# 论文图件（v2）：Figure 1 研究流程 / Figure 2 双队列森林图 / Figure 3 双队列一致性
# 输出：PDF（矢量）+ TIFF 300dpi + PNG 预览；英文标注
# ============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(cowplot) })
BASE <- getwd()
setwd(BASE)
OUT <- file.path(BASE, "figures"); dir.create(OUT, showWarnings = FALSE)
FONT <- "sans"

theme_pub <- theme_classic(base_size = 9, base_family = FONT) +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.35),
        axis.ticks = element_line(linewidth = 0.35),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold", size = 8.5),
        legend.title = element_blank(),
        plot.margin = margin(6, 10, 6, 6))

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".tiff")), p, width = w, height = h,
         dpi = 300, device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h,
         dpi = 300, bg = "white")
  cat("[OK]", name, "\n")
}

## ============================ Figure 1 ============================
boxes <- data.table(
  xmin = c(16,  2, 52,  8,  2, 51,  8, 18),
  xmax = c(84, 48, 98, 92, 49, 98, 92, 82),
  ymin = c(89, 72, 72, 55, 38, 38, 21,  5),
  ymax = c(96, 81, 81, 65, 48, 48, 31, 15),
  fill = c("#DCE9F5", "#E8F2E2", "#E8F2E2", "#F2F2F2", "#FDF2DC", "#FDF2DC", "#EDE7F6", "#F8DCDC"),
  size = c(2.7, 2.6, 2.6, 2.6, 2.7, 2.6, 2.6, 2.9),
  lab = c(
    "Exposure: gut microbiota (MiBioGen; n = 18,340)\n140 taxa \u00b7 1,748 instrument SNPs (1,363 unique)",
    "Discovery cohort: FinnGen R13 (DF13)\nAAA 4,634 \u00b7 TAA 5,704 \u00b7 AA 10,694 \u00b7 AD 1,219\n460,490 controls each",
    "Replication cohort: AAAgen (Nat Genet 2023)\n39,221 AAA cases / 1,086,107 controls\nhg19 coordinates; no rsIDs",
    "Harmonisation\nrsID matching (FinnGen) / hg19 position matching (AAAgen)\nallele alignment \u00b7 palindromic check \u00b7 1,343 AAAgen variants retained",
    "Two-sample MR\nIVW (primary) \u00b7 weighted median \u00b7 MR-Egger",
    "Sensitivity analyses\nCochran's Q \u00b7 MR-Egger intercept\nMR-PRESSO \u00b7 leave-one-out \u00b7 Steiger",
    "Replication criterion\nnominal P < 0.05 in both cohorts with concordant direction\nFDR q < 0.05 in both cohorts = robust association",
    "0 of 140 taxa replicated\nNo robust causal association")
)

arrows <- data.table(
  x  = c(30, 70, 25, 75, 30, 70, 25, 75, 50),
  y  = c(89, 89, 72, 72, 55, 55, 38, 38, 21),
  xe = c(30, 70, 25, 75, 30, 70, 25, 75, 50),
  ye = c(81, 81, 65, 65, 48, 48, 31, 31, 15)
)

p1 <- ggplot() +
  geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = boxes$fill, colour = "grey35", linewidth = 0.4) +
  geom_text(data = boxes, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
                              label = lab, size = size), lineheight = 1.18) +
  scale_size_identity() +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xe, yend = ye),
               arrow = arrow(length = unit(1.9, "mm"), type = "closed"),
               colour = "grey30", linewidth = 0.4) +
  annotate("text", x = 50, y = -2.8, size = 2.35, colour = "grey35", lineheight = 1.2,
           label = paste("Reverse MR (aortic disease \u2192 140 taxa; IVW) was performed for AA, AAA and TAA,",
                         "where \u2265 3 disease instruments were available (AD not estimable); all estimates were non-significant.",
                         sep = "\n")) +
  coord_cartesian(xlim = c(0, 100), ylim = c(-7, 98), expand = FALSE) +
  theme_void(base_size = 9, base_family = FONT) +
  theme(plot.margin = margin(4, 6, 4, 6))

save_fig(p1, "Figure1_study_design", 7.1, 6.9)

## ============================ Figure 2 ============================
xc <- fread(file.path(BASE, "results", "cross_cohort_lookup.csv"))
cand <- xc[p_R13 < 0.05 | p_AG < 0.05]
setorder(cand, p_R13)
cand[, taxon := factor(exposure, levels = rev(exposure))]
cat("[INFO] 森林图 taxa 数:", nrow(cand), "\n")

mk <- function(b, se, p, cohort) {
  data.table(taxon = cand$taxon, cohort = cohort,
             OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = p)
}
dd <- rbind(
  mk(cand$b_R13, cand$se_R13, cand$p_R13, "FinnGen R13\n(discovery)"),
  mk(cand$b_AG,  cand$se_AG,  cand$p_AG,  "AAAgen\n(replication)")
)
dd[, sig := p < 0.05]
dd[, lab := ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p))]
dd[, cohort := factor(cohort, levels = c("FinnGen R13\n(discovery)", "AAAgen\n(replication)"))]

p2 <- ggplot(dd, aes(x = OR, y = taxon)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.45, colour = "grey25") +
  geom_point(aes(fill = sig), shape = 21, size = 1.9, colour = "grey15", stroke = 0.55) +
  geom_text(aes(x = hi, label = lab), hjust = -0.30, vjust = 0.35, size = 2.3, colour = "grey20") +
  facet_wrap(~cohort, nrow = 1) +
  scale_x_log10(expand = expansion(mult = c(0.08, 0.32))) +
  scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "white")) +
  labs(x = "OR (95% CI) for abdominal aortic aneurysm\nper 1-SD higher taxon abundance", y = NULL) +
  theme_pub + theme(legend.position = "none", panel.spacing = unit(14, "pt"))

save_fig(p2, "Figure2_dual_cohort_forest", 7.1, 3.5)

## ============================ Figure 3 ============================
xc[, grp := fifelse(p_R13 < 0.05 | p_AG < 0.05, "nominal hit in \u2265 1 cohort", "neither")]
r <- cor(xc$b_R13, xc$b_AG)
tt <- r * sqrt((nrow(xc) - 2) / (1 - r^2))
pv <- 2 * pnorm(-abs(tt))

p3 <- ggplot(xc, aes(x = b_R13, y = b_AG)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey55", linewidth = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              linetype = "dashed", colour = "grey35", linewidth = 0.4) +
  geom_point(aes(colour = grp, shape = grp), size = 1.55, alpha = 0.9, stroke = 0.45) +
  scale_colour_manual(values = c("nominal hit in \u2265 1 cohort" = "#C0392B", "neither" = "grey35")) +
  scale_shape_manual(values = c("nominal hit in \u2265 1 cohort" = 17, "neither" = 16)) +
  labs(x = expression(beta ~ "(FinnGen R13, AAA)"), y = expression(beta ~ "(AAAgen, AAA)"),
       subtitle = sprintf("Pearson r = %.2f (P = %.3f); concordant direction %d/%d",
                          r, pv, xc[same_dir == TRUE, .N], nrow(xc))) +
  theme_pub +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(8, "pt"), legend.text = element_text(size = 7),
        plot.subtitle = element_text(size = 7.2, colour = "grey25"))

save_fig(p3, "Figure3_cross_cohort_concordance", 3.9, 3.9)

cat("\n========== 图件生成完成 ==========\n")
for (f in list.files(OUT)) cat(" ", f, "\n")
