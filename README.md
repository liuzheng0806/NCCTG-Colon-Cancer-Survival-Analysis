NCCTG 结肠癌辅助化疗试验 生存分析复盘项目。

## 依赖

- R(版本 4.0 以上,本项目实测于 4.3.3)
- R 包:`survival`(base R 安装时通常自带,如无则 `install.packages("survival")`)

**仅此一项依赖。** 不依赖 dplyr、tidyverse、smcfcs 等额外包。

## 运行顺序

六个脚本严格按顺序运行。后一个脚本依赖前一个产出的 `adtte.RData`。

| 顺序 | 脚本 | 输入 | 输出 |
|---|---|---|---|
| 1 | `3A_data_setup.R` | (自带 `survival::colon`) | `adtte.RData` |
| 2 | `3B_table1.R` | `adtte.RData` | 控制台 Table 1 |
| 3 | `3C_primary_analysis.R` | `adtte.RData` | 控制台主分析结果 + 图形窗口 KM 曲线 |
| 4 | `3D_ph_diagnostics.R` | `adtte.RData` | 控制台 cox.zph + 图形窗口 Schoenfeld 残差图 |
| 5 | `3E_sensitivity.R` | `adtte.RData` | 控制台 SA-1/SA-2/SA-3 结果 + 图形窗口 tipping point |
| 6 | `3F_subgroup_forest.R` | `adtte.RData` | 控制台亚组表 + 图形窗口 forest plot |

## 运行方式

在 R 或 RStudio 中,把工作目录设到脚本所在文件夹,然后:

```r
setwd("脚本所在目录")
source("3A_data_setup.R")
source("3B_table1.R")
source("3C_primary_analysis.R")
source("3D_ph_diagnostics.R")
source("3E_sensitivity.R")
source("3F_subgroup_forest.R")
```

绘图脚本(3C / 3D / 3E / 3F)会直接在 R/RStudio 图形窗口绘图,不写文件。如需存 PDF,在脚本绘图段前后加 `pdf()` / `dev.off()`。

## 可重复性

- 涉及随机抽样的脚本(`3E_sensitivity.R` 中 SA-3 部分)已固定随机种子(`set.seed(20250526)`),多次运行结果一致。
- 主分析所有方法参数(Cox 的 tied events 用 Efron、KM 的 CI 用 log-log 等)在脚本中显式声明,不依赖 R 默认值。

## 数据来源

`survival::colon`,R `survival` 包自带。原始研究为 Moertel et al., NEJM 1990, 322:352-358。本项目使用 ITT 人群 N=619(Obs 与 Lev+5FU 两臂,剔除 Lev 单药),终点为无复发生存期(RFS, etype=1)。
