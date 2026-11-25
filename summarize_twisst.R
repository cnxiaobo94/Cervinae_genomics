
library(data.table)

# 1. 设置文件路径和前缀
setwd("D:/Workspace/Cervinae/wg_analyses/5_phytop/twisst")
prefix <- "Cervinae_10inds_1out_auto_w1000kb_S1000kb_phyML"
outprefix <- "twisst_ElaDav_CerAlb_CerNip_Reddeer_w1000kb_S1000kb_phyML"

# Cervinae_10inds_1out_auto_w50kb_S50kb_phyML
# Cervinae_10inds_1out_auto_w100kb_S100kb_phyML
# Cervinae_10inds_1out_auto_w500kb_S500kb_phyML
# Cervinae_10inds_1out_auto_w1000kb_S1000kb_phyML

# twisst_CerAlb_CerEla_CerHan_CerCanAS
# twisst_CerNip_CerAlb_CerEla_CerHanCerCanAS
# twisst_DamDam_ElaDav_RusAlfRucEld_Cervus

data_file <- paste0(prefix, ".data.tsv")
tree_file <- paste0(prefix, ".trees.gz")
weight_file <- paste0(outprefix, ".weights.csv")
outdir <- paste0(outprefix, "_twisst_summary")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

topos <- c("topo1", "topo2", "topo3")

# 2. 读取窗口坐标
d <- fread(data_file, na.strings = c("NA", "nan", "NaN", ""))
required <- c("scaffold", "start", "end", "sites", "lnL")
if (!all(required %in% names(d))) stop("data.tsv 缺少必要列")

d[, `:=`(
  start = as.numeric(start),
  end = as.numeric(end),
  sites = as.numeric(sites),
  lnL = suppressWarnings(as.numeric(lnL)),
  window_id = .I
)]

# 3. 读取树文件并识别有效 Newick 树
con <- gzfile(tree_file, "rt")
trees <- readLines(con, warn = FALSE)
close(con)

is_tree <- function(x) {
  !is.na(x) &&
    grepl("\\(", x) &&
    grepl(";", x) &&
    !grepl("^\\s*NA\\s*$", x, ignore.case = TRUE)
}

tree_ok <- vapply(trees, is_tree, logical(1))

if (length(trees) == nrow(d)) {
  tree_ok_window <- tree_ok
} else if (all(tree_ok) && length(trees) == sum(!is.na(d$lnL))) {
  tree_ok_window <- rep(FALSE, nrow(d))
  tree_ok_window[which(!is.na(d$lnL))] <- TRUE
} else {
  stop("树文件与 data.tsv 的窗口无法可靠匹配，请检查 NA 行和窗口记录")
}

# 4. 读取 Twisst 权重
read_weights <- function(path) {
  x <- fread(path, header = TRUE)
  if (ncol(x) < 3) stop("权重文件少于三列")
  
  w <- as.data.frame(lapply(x[, 1:3], as.numeric))
  names(w) <- topos
  w
}

w <- read_weights(weight_file)

# 5. 将权重对应到窗口
if (nrow(w) == nrow(d)) {
  weight_window_id <- seq_len(nrow(d))
} else if (nrow(w) == sum(tree_ok_window)) {
  weight_window_id <- which(tree_ok_window)
} else {
  stop("Twisst 权重行数无法与窗口对应，请检查权重文件是否含表头或 NA 窗口")
}

d[, `:=`(
  tree_valid = tree_ok_window,
  topo1 = NA_real_,
  topo2 = NA_real_,
  topo3 = NA_real_
)]

d[weight_window_id, `:=`(
  topo1 = w[[1]],
  topo2 = w[[2]],
  topo3 = w[[3]]
)]

# 6. 计算每个窗口的权重比例
weight_mat <- as.matrix(d[, ..topos])
weight_sum <- rowSums(weight_mat)

d[, weight_valid := tree_valid &
    apply(weight_mat, 1, function(x) all(is.finite(x))) &
    is.finite(weight_sum) & weight_sum > 0]

for (t in topos) {
  prop_col <- paste0(t, "_prop")
  d[, (prop_col) := NA_real_]
  d[weight_valid == TRUE,
    (prop_col) := get(t) / rowSums(.SD),
    .SDcols = topos]
}

# 7. 确定 winner，权重并列不指定 winner
d[, `:=`(winner = NA_character_, winner_tie = FALSE)]

for (i in which(d$weight_valid)) {
  z <- as.numeric(d[i, ..topos])
  winners <- which(abs(z - max(z)) < 1e-10)
  
  if (length(winners) == 1) {
    d$winner[i] <- topos[winners]
  } else {
    d$winner_tie[i] <- TRUE
  }
}

fwrite(d, file.path(outdir, "window_results.tsv"), sep = "\t")

# 8. 全基因组拓扑权重汇总，增加 SD 和四分位数
prop_cols <- paste0(topos, "_prop")
valid <- d$weight_valid

topology_summary <- rbindlist(lapply(seq_along(topos), function(i) {
  x <- d[[prop_cols[i]]][valid]
  x <- x[is.finite(x)]
  
  data.table(
    topology = topos[i],
    n_valid_windows = length(x),
    mean_weight_proportion = if (length(x)) mean(x) else NA_real_,
    sd_weight_proportion = if (length(x) > 1) sd(x) else NA_real_,
    median_weight_proportion = if (length(x)) median(x) else NA_real_,
    q25_weight_proportion = if (length(x)) as.numeric(quantile(x, 0.25)) else NA_real_,
    q75_weight_proportion = if (length(x)) as.numeric(quantile(x, 0.75)) else NA_real_,
    winner_windows = sum(d$winner == topos[i], na.rm = TRUE),
    winner_fraction = if (sum(valid) > 0) {
      sum(d$winner == topos[i], na.rm = TRUE) / sum(valid)
    } else {
      NA_real_
    }
  )
}))

fwrite(topology_summary,
       file.path(outdir, "topology_summary.tsv"), sep = "\t")

# 9. 基础统计
stats <- data.table(
  metric = c("total_windows", "valid_tree_windows",
             "valid_weight_windows", "invalid_tree_windows",
             "invalid_weight_windows", "winner_tie_windows",
             "winner_defined_windows"),
  value = c(nrow(d), sum(d$tree_valid, na.rm = TRUE),
            sum(valid, na.rm = TRUE),
            sum(!d$tree_valid, na.rm = TRUE),
            sum(d$tree_valid & !valid, na.rm = TRUE),
            sum(d$winner_tie, na.rm = TRUE),
            sum(!is.na(d$winner)))
)

fwrite(stats, file.path(outdir, "basic_statistics.tsv"), sep = "\t")

# 10. 统计连续 winner 区段
# 无效窗口、winner 并列、scaffold 变化、坐标不连续均中断区段
setorder(d, scaffold, start, end)

runs <- list()
current <- integer(0)

save_run <- function(idx) {
  if (!length(idx)) return(NULL)
  
  x <- d[idx]
  data.table(
    scaffold = x$scaffold[1],
    topology = x$winner[1],
    start = min(x$start),
    end = max(x$end),
    windows = nrow(x),
    length_bp = max(x$end) - min(x$start) + 1
  )
}

for (i in seq_len(nrow(d))) {
  usable <- isTRUE(d$weight_valid[i]) && !is.na(d$winner[i])
  
  if (!usable) {
    if (length(current)) {
      runs[[length(runs) + 1]] <- save_run(current)
      current <- integer(0)
    }
    next
  }
  
  if (!length(current)) {
    current <- i
    next
  }
  
  j <- tail(current, 1)
  
  adjacent <- identical(as.character(d$scaffold[i]),
                        as.character(d$scaffold[j])) &&
    identical(d$winner[i], d$winner[j]) &&
    d$start[i] == d$end[j] + 1
  
  if (adjacent) {
    current <- c(current, i)
  } else {
    runs[[length(runs) + 1]] <- save_run(current)
    current <- i
  }
}

if (length(current)) {
  runs[[length(runs) + 1]] <- save_run(current)
}

run_df <- rbindlist(runs, fill = TRUE)

if (!nrow(run_df)) {
  run_df <- data.table(
    scaffold = character(),
    topology = character(),
    start = numeric(),
    end = numeric(),
    windows = integer(),
    length_bp = numeric()
  )
}

fwrite(run_df, file.path(outdir, "topology_runs.tsv"), sep = "\t")

# 11. 连续区段长度汇总，增加 SD 和四分位数
run_summary <- rbindlist(lapply(topos, function(t) {
  x <- run_df[topology == t, length_bp]
  x <- x[is.finite(x)]
  
  data.table(
    topology = t,
    winner_windows = sum(d$winner == t, na.rm = TRUE),
    run_count = length(x),
    mean_length_bp = if (length(x)) mean(x) else NA_real_,
    sd_length_bp = if (length(x) > 1) sd(x) else NA_real_,
    median_length_bp = if (length(x)) median(x) else NA_real_,
    q25_length_bp = if (length(x)) as.numeric(quantile(x, 0.25)) else NA_real_,
    q75_length_bp = if (length(x)) as.numeric(quantile(x, 0.75)) else NA_real_,
    max_length_bp = if (length(x)) max(x) else NA_real_
  )
}))

fwrite(run_summary, file.path(outdir, "run_summary.tsv"), sep = "\t")

message("分析完成，结果保存在：", outdir)
