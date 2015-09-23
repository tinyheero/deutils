#' Run a Wilcox Test
#'
#' @param vals Numeric vector with sample names as keys
#' @param groupA Sample names from Group A
#' @param groupB Sample names from Group B
#' @param pairedFlag Run a paired wilcox mode
#' @export
run_de_test <- function(vals, groupA, groupB, pairedFlag = FALSE) {
  res <- try(wilcox.test(x = vals[groupA], y = vals[groupB], 
                         paired = pairedFlag))

  if (pairedFlag){
    fc.diff <- mean(vals[groupA] - vals[groupB])
  } else{
    fc.diff <- mean(vals[groupA]) - mean(vals[groupB])
  }

  outList <- list("meanExprsA" = mean(vals[groupA]), 
                  "meanExprsB" = mean(vals[groupB]), 
                  "fcDiff" = fc.diff)

  if (class(res) != "try-error"){
    outList[["pval"]] <- res$p.value
  } else{
    outList[["pval"]] <- NA
  }

  outList
}
