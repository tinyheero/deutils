#' Run a Differential Expression Test
#'
#' @param gene.exprs.mat Gene expression matrix with samples as columns and
#'  genes as rows
#' @param sample.annot.df Sample annotation file that must contain the columns
#'  sampleID and condition. The condition column should be a factor with the 
#'  levels determining which groups. It will be level 1 - level 2
#' @param pairedFlag Run a paired wilcox mode
#' @export
run_de <- function(gene.exprs.mat, sample.annot.df, test.type = "wilcox",
                   log.flag = TRUE) {
  condition1 <- levels(sample.annot.df[["condition"]])[1]
  condition2 <- levels(sample.annot.df[["condition"]])[2]

  message(paste0("Differential Expression: ", condition1, " vs. ", 
                 condition2))

  filter_crit <- lazyeval::interp(~ condition == condition1, 
                                  condition = as.name("condition"))
  condition1.samples <- dplyr::filter_(sample.annot.df, 
                                       filter_crit)[["sampleID"]]

  filter_crit <- lazyeval::interp(~ condition == condition2, 
                                  condition = as.name("condition"))
  condition2.samples <- dplyr::filter_(sample.annot.df, 
                                       filter_crit)[["sampleID"]]

  samples <- c(condition1.samples, condition2.samples)
  wilcox.res <- apply(gene.exprs.mat, 1, run_wilcox_test, condition1.samples, 
                      condition2.samples)

  wilcox.res.df <- wilcox.res %>%
    lapply(dplyr::as_data_frame) %>% 
    dplyr::bind_rows(.id = "geneID") %>%
    dplyr::rename_(.dots = setNames(list("meanExprsA", "meanExprsB"),
                                    c(paste0("meanExprs", condition1),
                                      paste0("meanExprs", condition2))))

  wilcox.res.df
}
