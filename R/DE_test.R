library(tidyverse)

f_ls <- list.files("~/Documents/Exported Data c", full.names = T)

f_name <- sapply(strsplit(f_ls, split = "/"), function(x) x[6])

f_pop <- sapply(strsplit(f_name, split = "_"), function(x){
  temp <- gsub("_NA", "", paste(x[6], x[7], sep = "_"))
  gsub(".csv", "", temp)
})


f_sample <- sapply(strsplit(f_name, split = "_"), function(x) x[2])
f_timepoint <- sapply(strsplit(f_name, split = "_"), function(x) x[3])
f_rep <- sapply(strsplit(f_name, split = "_"), function(x) x[4])

col_data <- data.frame(sample = f_sample,
                       timepoint = f_timepoint,
                       #rep = f_rep,
                       #population = f_pop,
                       rep_pop = paste(f_pop, f_rep, sep = ";"))

csv_ls <- lapply(f_ls, function(x) read.csv(x))
summary(csv_ls[[1]])

colnames(csv_ls[[1]])
marker_stat <- c("Ki67", "PD.L1", "PD1", "GzB", "CD25", "CD38")

#median all
med <- t(sapply(csv_ls, function(x) (apply(x, 2, median))))

#cellcount
pop_counts <- cbind(col_data,
                    counts = sapply(csv_ls, function(x) nrow(x))) %>%
  spread(key = "rep_pop", value = "counts")

n_pop <- length(unique(col_data$rep_pop))
design <- pop_counts[, 1 : (ncol(pop_counts)- n_pop)]
pop_counts <- t(pop_counts[, (ncol(pop_counts)- n_pop + 1) : ncol(pop_counts)])

med_ls <- list()
for (i in marker_stat) {
  temp <- cbind(col_data, 
        marker = med[, i]) %>%
    spread(key = "rep_pop", value = "marker") 
  n_pop <- length(unique(col_data$rep_pop))
  med_ls[[i]] <- t(temp[, (ncol(temp)- n_pop + 1) : ncol(temp)])
}

colnames(design[[1]])
formula <- as.formula(paste("~", paste(c("sample", "timepoint"), collapse = " + ")))
design <- model.matrix(formula, data = design)

contrast <- matrix(c(0, 0, 0, 1), ncol = 1)

library(limma)

res <- list()
for (i in names(med_ls)) {
  dupcor <- duplicateCorrelation(med_ls[[1]], design = design, ndups = 2, spacing = 1)
  
  fit <- lmFit(med_ls[[1]], design, weights = pop_counts, ndups = 2,correlation = dupcor$consensus.correlation)
  
  fit <- contrasts.fit(fit, contrast)
  
  efit <- eBayes(fit)
  
  res[[i]] <- topTable(efit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "none")
  rownames(res[[i]]) <- unique(f_pop)
}


library(edgeR)
norm_factors <- calcNormFactors(pop_counts, method = "TMM")
y <- DGEList(pop_counts, norm.factors = norm_factors)
y <- estimateDisp(y, design, trend.method = NULL)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, contrast = contrast)
top <- edgeR::topTags(lrt, n = Inf, adjust.method = "BH", sort.by = "none")
