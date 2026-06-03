r# ============================================================
# 项目：NCCTG结肠癌辅助化疗试验 生存分析复盘
# 脚本：3-A 环境搭建 + 数据加载 + ITT人群构造
# 数据来源：survival::colon (R package)
# 对应SAP章节：§3 Study Design, §4 Analysis Populations
# ============================================================
# ---- 0. 环境准备 --------------------------------------------
# 清空工作环境，保证可重复性
rm(list = ls())
# 记录R版本（对应SAP §9 Software & Reproducibility）
# 在真实SAP里这个版本号要写死
cat("R version:\n")
print(R.version.string)
# 加载必要包
# 只依赖 survival（base R），避免额外依赖
library(survival)
cat("\nsurvival 包版本:", as.character(packageVersion("survival")), "\n")
# ---- 1. 加载原始数据 ----------------------------------------
# colon数据集说明：
#   - 每个患者有2行，由etype区分：
#     etype = 1: recurrence（复发事件行）
#     etype = 2: death（死亡事件行）
#   - 真 RFS 终点 = 复发 或 全因死亡，先发生者；需合并这两行构造（见第1节末）
#   - 治疗组rx：Obs / Lev / Lev+5FU（三臂）
data(cancer, package = "survival")  # colon在cancer里
# 核查原始数据结构
cat("\n=== 原始数据结构 ===\n")
cat("总行数（应为929患者×2行 = 1858）:", nrow(colon), "\n")
cat("变量列表:\n")
print(names(colon))
cat("\n治疗组×etype 分布（原始三臂）:\n")
print(table(colon$rx, colon$etype,
            dnn = c("治疗组", "etype(1=recurrence,2=death)")))
# ===== 真 RFS 终点构造：合并 etype 两行 =====
# 关键：RFS = recurrence OR all-cause death，先发生者。
#   - 只取 etype==1 单行 会把"无复发即死亡"错当删失 → 得到 time-to-recurrence，
#     而非 RFS。本项目曾踩此坑并修正，这里显式合并两行。
rec   <- colon[colon$etype == 1, ]   # 每患者的复发行
death <- colon[colon$etype == 2, ]   # 每患者的死亡行
rec   <- rec[order(rec$id), ]
death <- death[order(death$id), ]
stopifnot(all(rec$id == death$id))   # 两行 id 必须一一对应
rfs <- rec  # 以复发行为骨架（协变量、rx、id 都在这一行）
# 事件指示：复发(rec$status==1) 或 死亡(death$status==1) 任一发生即为事件
rfs$status <- as.integer(rec$status == 1 | death$status == 1)
# 事件/删失时间：
#   复发者 → 取复发时间 rec$time
#   未复发者 → 取 death 行时间（死亡时间，或无事件时的随访删失时间）
rfs$time <- ifelse(rec$status == 1, rec$time, death$time)
# 核对（全三臂）：RFS 事件应为 Obs 190 / Lev 182 / Lev+5FU 134
cat("三臂真RFS事件:\n"); print(table(rfs$rx, rfs$status))
# ---- 2. 构造ITT分析数据集 -----------------------------------
# 对应SAP §4：ITT Population
# 规则：
#   (a) 数据来自上面合并好的真 RFS（rfs），不是原始单行
#   (b) 只保留 Obs 和 Lev+5FU 两个臂（剔除 Lev 单药组）
#   (c) 按SAP预先指定，Obs 为参照组
# (b) 在合并后的 rfs 上筛两臂
adtte <- rfs[rfs$rx %in% c("Obs", "Lev+5FU"), ]
# (c) 显式设定因子水平和参照组
#     第一个水平 = Cox参照组；Obs在前 → HR解读为"Lev+5FU相对于Obs的风险比"
adtte$rx <- factor(adtte$rx, levels = c("Obs", "Lev+5FU"))

# 派生变量1：node4（是否>4个阳性淋巴结）
# 用途：SAP预设的分层因子
adtte$node4 <- factor(
  ifelse(adtte$nodes > 4, ">4 nodes", "<=4 nodes"),
  levels = c("<=4 nodes", ">4 nodes")
)

# 派生变量2：时间单位转换 天→月
# 用途：报中位生存时间时用月更直观（临床习惯）
# 转换系数：365.25/12 ≈ 30.4375天/月
adtte$time_months <- adtte$time / (365.25 / 12)

# 派生变量3：differ转因子（肿瘤分化程度）
# 原始编码：1=well, 2=moderate, 3=poor
adtte$differ <- factor(adtte$differ,
                       levels = c(1, 2, 3),
                       labels = c("Well", "Moderate", "Poor"))

# 派生变量4：extent转因子（肿瘤局部侵犯程度）
# 原始编码：1=submucosa, 2=muscle, 3=serosa, 4=contiguous
adtte$extent <- factor(adtte$extent,
                       levels = c(1, 2, 3, 4),
                       labels = c("Submucosa", "Muscle",
                                  "Serosa", "Contiguous"))

# 派生变量5：sex转因子
adtte$sex <- factor(adtte$sex,
                    levels = c(0, 1),
                    labels = c("Female", "Male"))


# ---- 3. 数据核查（Data Verification）-----------------------
# 对应SAP §9：Reproducibility
# 这一步必须跑，确认人群构造结果符合预期

cat("\n=== ITT人群核查 ===\n")
cat("ITT总N（预期619）:", nrow(adtte), "\n")

cat("\n治疗组分布:\n")
print(table(adtte$rx, dnn = "治疗组"))

cat("\n事件/删失分布:\n")
print(table(adtte$rx, adtte$status,
            dnn = c("治疗组", "status(1=事件,0=删失)")))

cat("\n分层因子node4分布:\n")
print(table(adtte$rx, adtte$node4,
            dnn = c("治疗组", "node4")))

cat("\n=== 关键协变量缺失情况 ===\n")
miss_vars <- c("nodes", "differ", "extent")
for (v in miss_vars) {
  n_miss <- sum(is.na(adtte[[v]]))
  cat(sprintf("  %-8s: NA = %d 人 (%.1f%%)\n",
              v, n_miss, 100 * n_miss / nrow(adtte)))
}

# 缺失核查（仅描述，不构造独立人群）
# 注意：主分析模型 rx+strata(node4)+age+sex+extent 只用到 nodes/extent/age/sex，
#       其中仅 nodes 含缺失(12)，故主分析实际剔除 12 人、拟合 N=607。
#       differ(缺13) 不在主模型，不影响主分析 N。
cat("\n各协变量缺失：nodes=", sum(is.na(adtte$nodes)),
    " differ=", sum(is.na(adtte$differ)),
    " extent=", sum(is.na(adtte$extent)),
    "（differ 不入主模型，不影响主分析 N）\n", sep = "")
cat("主分析实际拟合 N（仅因 nodes 缺失剔 12）:", 
    sum(!is.na(adtte$nodes)), "\n")




# ---- 4. 保存分析数据集 --------------------------------------
# 在真实项目里，分析数据集要锁定（lock）后才能使用
# 这里我们保存为.RData，后续脚本直接load

save(adtte, file = "adtte.RData")
cat("\n数据集已保存：adtte.RData\n")
cat("  adtte: ITT人群（真RFS），N =", nrow(adtte), "\n")

cat("\n=== 3-A 完成 ===\n")
