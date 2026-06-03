# ============================================================
# 项目：NCCTG结肠癌辅助化疗试验 生存分析复盘
# 脚本：3-B Table 1 基线特征表
# 对应SAP章节：§5 报告内容 / TFL: Table 1
# 输入：adtte.RData（由3-A生成）
# ============================================================

rm(list = ls())
library(survival)

# ---- 0. 载入3-A构造的分析数据集 -----------------------------
load("adtte.RData")
cat("ITT人群 N =", nrow(adtte), "\n\n")


# ---- 1. 辅助函数：格式化输出 --------------------------------

# 连续变量：返回 "mean ± SD" 字符串
fmt_cont <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.1f ± %.1f", mean(x), sd(x))
}

# 分类变量某一水平：返回 "n (xx.x%)" 字符串
# 百分比分母用非缺失样本数（这是标准做法）
fmt_cat <- function(x, level) {
  x_nonmiss <- x[!is.na(x)]
  n <- sum(x_nonmiss == level)
  pct <- 100 * n / length(x_nonmiss)
  sprintf("%d (%.1f%%)", n, pct)
}

# 缺失数：返回 "n" 或在无缺失时返回 "0"
fmt_miss <- function(x) {
  as.character(sum(is.na(x)))
}


# ---- 2. 按治疗臂拆分数据 ------------------------------------
obs   <- adtte[adtte$rx == "Obs", ]
combo <- adtte[adtte$rx == "Lev+5FU", ]
all   <- adtte


# ---- 3. 逐行构造Table 1 -------------------------------------
# 每一行是一个 c(变量名, Obs列, Lev+5FU列, Overall列)

cat("============================================================\n")
cat("  Table 1. Baseline Characteristics by Treatment Arm\n")
cat("  (ITT Population, N = 619)\n")
cat("============================================================\n")

# 表头
hdr <- sprintf("%-32s %-18s %-18s %-18s",
               "Characteristic",
               sprintf("Obs (N=%d)", nrow(obs)),
               sprintf("Lev+5FU (N=%d)", nrow(combo)),
               sprintf("Overall (N=%d)", nrow(all)))
cat(hdr, "\n")
cat(strrep("-", 88), "\n")

# 行输出函数
prow <- function(label, c_obs, c_combo, c_all) {
  cat(sprintf("%-32s %-18s %-18s %-18s\n",
              label, c_obs, c_combo, c_all))
}

# --- Age（连续变量）---
prow("Age (years), mean ± SD",
     fmt_cont(obs$age), fmt_cont(combo$age), fmt_cont(all$age))

# --- Sex（分类）---
cat("Sex, n (%)\n")
prow("  Male",
     fmt_cat(obs$sex,"Male"), fmt_cat(combo$sex,"Male"),
     fmt_cat(all$sex,"Male"))
prow("  Female",
     fmt_cat(obs$sex,"Female"), fmt_cat(combo$sex,"Female"),
     fmt_cat(all$sex,"Female"))

# --- node4（分层因子）---
cat("Positive nodes, n (%)\n")
prow("  <=4 nodes",
     fmt_cat(obs$node4,"<=4 nodes"), fmt_cat(combo$node4,"<=4 nodes"),
     fmt_cat(all$node4,"<=4 nodes"))
prow("  >4 nodes",
     fmt_cat(obs$node4,">4 nodes"), fmt_cat(combo$node4,">4 nodes"),
     fmt_cat(all$node4,">4 nodes"))

# --- differ（肿瘤分化，有缺失）---
cat("Differentiation, n (%)\n")
for (lv in c("Well","Moderate","Poor")) {
  prow(paste0("  ", lv),
       fmt_cat(obs$differ,lv), fmt_cat(combo$differ,lv),
       fmt_cat(all$differ,lv))
}
prow("  Missing",
     fmt_miss(obs$differ), fmt_miss(combo$differ),
     fmt_miss(all$differ))

# --- extent（肿瘤侵犯程度）---
cat("Extent of invasion, n (%)\n")
for (lv in c("Submucosa","Muscle","Serosa","Contiguous")) {
  prow(paste0("  ", lv),
       fmt_cat(obs$extent,lv), fmt_cat(combo$extent,lv),
       fmt_cat(all$extent,lv))
}

# --- obstruct（结肠梗阻，0/1变量，作为附加基线特征）---
cat("Colon obstruction, n (%)\n")
prow("  Yes",
     fmt_cat(obs$obstruct,1), fmt_cat(combo$obstruct,1),
     fmt_cat(all$obstruct,1))

# --- 事件数（不是基线特征，但放表末供参考）---
cat(strrep("-", 88), "\n")
prow("Recurrence/death events, n (%)",
     fmt_cat(obs$status,1), fmt_cat(combo$status,1),
     fmt_cat(all$status,1))

cat("============================================================\n")
cat("注：百分比分母为非缺失样本数。\n")
cat("    按CONSORT规范，随机化试验Table 1不报组间检验p值。\n")
cat("============================================================\n")

cat("\n=== 3-B 完成 ===\n")
