# ============================================================
# 项目：NCCTG结肠癌辅助化疗试验 生存分析复盘
# 脚本：3-F 亚组分析 + forest plot
# 对应SAP章节：§5.4 Subgroup Analyses
# 输入：adtte.RData（由3-A生成）
#
# ⚠️ 设计选择：
#   - 各亚组内 HR：仅 rx 的 Cox（不调整协变量），
#     看该亚组粗治疗效应，保持 forest plot 干净
#   - 交互作用 p 值：用含 rx×亚组 交互项的全 ITT 模型
#   - 所有亚组分析均为 exploratory，不做 confirmatory 主张
#   - 亚组内不报 p 值，只报 HR + 95%CI（多重比较考量）
#   本脚本所有部分均为 base R，已实测验证。
# ============================================================

rm(list = ls())
library(survival)
load("adtte.RData")


# ============================================================
# 第一部分：定义预设亚组
# ============================================================
# 派生亚组变量（部分3-A已有，这里统一构造二分类）

adtte$sg_age <- factor(ifelse(adtte$age < 60, "Age <60", "Age >=60"),
                       levels = c("Age <60", "Age >=60"))
adtte$sg_sex <- adtte$sex                          # 已是factor
adtte$sg_node <- adtte$node4                        # 已是 <=4 / >4
# differ 二分类：Well/Moderate vs Poor
adtte$sg_diff <- factor(
  ifelse(adtte$differ == "Poor", "Poorly diff.", "Well/Mod diff."),
  levels = c("Well/Mod diff.", "Poorly diff."))
# extent 二分类：<=Muscle vs >Muscle(Serosa/Contiguous)
adtte$sg_ext <- factor(
  ifelse(adtte$extent %in% c("Submucosa", "Muscle"),
         "<=Muscle", ">Muscle"),
  levels = c("<=Muscle", ">Muscle"))

# 亚组变量清单：显示名 -> 列名
subgroups <- list(
  "Age"           = "sg_age",
  "Sex"           = "sg_sex",
  "Positive nodes"= "sg_node",
  "Differentiation"="sg_diff",
  "Extent"        = "sg_ext"
)


# ============================================================
# 第二部分：计算各亚组内的治疗 HR
# ============================================================
# 每个亚组的每个水平，单独跑 Cox(仅rx)，提取HR和95%CI

# 存放forest plot数据：label, n, events, HR, lower, upper
fp <- data.frame(label = character(), n = integer(),
                 events = integer(), HR = numeric(),
                 lower = numeric(), upper = numeric(),
                 is_header = logical(),
                 stringsAsFactors = FALSE)

# 先放Overall(全人群)
fit_all <- coxph(Surv(time_months, status) ~ rx,
                 data = adtte, ties = "efron")
ci_all <- summary(fit_all)$conf.int["rxLev+5FU", ]
fp <- rbind(fp, data.frame(
  label = "Overall", n = nrow(adtte),
  events = sum(adtte$status),
  HR = ci_all["exp(coef)"], lower = ci_all["lower .95"],
  upper = ci_all["upper .95"], is_header = FALSE))

# 逐个亚组
for (sg_name in names(subgroups)) {
  col <- subgroups[[sg_name]]
  
  # 亚组标题行（无HR）
  fp <- rbind(fp, data.frame(
    label = sg_name, n = NA, events = NA,
    HR = NA, lower = NA, upper = NA, is_header = TRUE))
  
  # 该亚组各水平
  for (lv in levels(adtte[[col]])) {
    sub <- adtte[!is.na(adtte[[col]]) & adtte[[col]] == lv, ]
    # 该水平内单独跑 Cox(仅rx)
    fit <- tryCatch(
      coxph(Surv(time_months, status) ~ rx, data = sub, ties = "efron"),
      error = function(e) NULL)
    if (!is.null(fit)) {
      ci <- summary(fit)$conf.int["rxLev+5FU", ]
      fp <- rbind(fp, data.frame(
        label = paste0("  ", lv), n = nrow(sub),
        events = sum(sub$status),
        HR = ci["exp(coef)"], lower = ci["lower .95"],
        upper = ci["upper .95"], is_header = FALSE))
    }
  }
}
rownames(fp) <- NULL


# ============================================================
# 第三部分：交互作用检验
# ============================================================
# 用含 rx×亚组 交互项的模型，似然比检验交互项是否显著
# 这才是"治疗效应是否随亚组变化"的正确问法

cat("=== 交互作用检验（治疗 × 亚组）===\n")
cat("(评估治疗效应是否在亚组间异质；非各亚组分别看p值)\n\n")

interaction_p <- list()
for (sg_name in names(subgroups)) {
  col <- subgroups[[sg_name]]
  dat <- adtte[!is.na(adtte[[col]]), ]
  
  # 无交互模型 vs 有交互模型，似然比检验
  f0 <- coxph(Surv(time_months, status) ~ rx + dat[[col]],
              data = dat, ties = "efron")
  f1 <- coxph(Surv(time_months, status) ~ rx * dat[[col]],
              data = dat, ties = "efron")
  lrt <- anova(f0, f1)
  p_int <- lrt$`Pr(>|Chi|)`[2]
  interaction_p[[sg_name]] <- p_int
  
  cat(sprintf("  %-18s 交互 p = %.3f  %s\n",
              sg_name, p_int,
              ifelse(p_int < 0.05,
                     "<- 提示效应可能异质(探索性)", "")))
}


# ============================================================
# 第四部分：打印 forest plot 数据表
# ============================================================
cat("\n=== Forest Plot 数据（Table: 亚组分析）===\n")
cat(sprintf("%-22s %6s %7s   %s\n",
            "Subgroup", "N", "Events", "HR (95% CI)"))
cat(strrep("-", 60), "\n")
for (i in seq_len(nrow(fp))) {
  if (fp$is_header[i]) {
    cat(sprintf("%-22s\n", fp$label[i]))
  } else {
    cat(sprintf("%-22s %6s %7s   %.2f (%.2f-%.2f)\n",
                fp$label[i],
                ifelse(is.na(fp$n[i]), "", fp$n[i]),
                ifelse(is.na(fp$events[i]), "", fp$events[i]),
                fp$HR[i], fp$lower[i], fp$upper[i]))
  }
}


# ============================================================
# 第五部分：绘制 forest plot（窗口输出）
# ============================================================
# 横轴：HR(对数刻度)；每行一个亚组水平：方块=HR，横线=95%CI
# 竖虚线 HR=1（无效应线）

plot_rows <- fp[!fp$is_header, ]   # 只画有HR的行
n_all <- nrow(fp)

# 画布：左侧留标签空间
par(mar = c(4.5, 11, 3, 2))
ylim <- c(1, n_all)
xlim <- c(0.2, 2.0)   # HR范围(对数轴)

plot(NA, xlim = xlim, ylim = c(n_all, 1), log = "x",
     xlab = "Hazard Ratio (Lev+5FU vs Obs)", ylab = "",
     yaxt = "n", main = "Figure 2. Subgroup Analysis (Forest Plot)")

# HR=1 无效应线
abline(v = 1, lty = 2, col = "gray40")

# 逐行画
for (i in seq_len(nrow(fp))) {
  y <- i
  if (fp$is_header[i]) {
    # 亚组标题：左对齐加粗显示
    mtext(fp$label[i], side = 2, at = y, las = 1,
          line = 10, adj = 0, font = 2, cex = 0.8)
  } else {
    # 标签
    mtext(fp$label[i], side = 2, at = y, las = 1,
          line = 10, adj = 0, cex = 0.75)
    # CI 横线
    lines(c(fp$lower[i], fp$upper[i]), c(y, y), col = "steelblue", lwd = 1.5)
    # HR 方块（Overall用实心大方块）
    pch_i <- if (fp$label[i] == "Overall") 18 else 15
    cex_i <- if (fp$label[i] == "Overall") 1.6 else 1.1
    points(fp$HR[i], y, pch = pch_i, col = "firebrick", cex = cex_i)
  }
}

cat("\nForest plot 已绘制到图形窗口。\n")
cat("判读：方块=HR点估计，横线=95%CI，竖虚线=HR=1(无效应)。\n")
cat("      CI 跨过1的亚组，治疗效应在该组不确定。\n")
cat("      看效应方向是否一致，不挑单个'显著'亚组。\n")

cat("\n=== 3-F 完成 ===\n")
