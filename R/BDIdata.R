#' BDIdata dataset
#'
#' This dataset includes 557 depressed patients (total 7117 observations) in the cognitive
#' behavior therapy arm in the Enhancing Recovery in Coronary Heart Disease Patients (ENRICHD)
#' study.
#'
#' \itemize{
#'   \item ID.        Subject ID
#'   \item time.      Study visit time (in days) since randomization
#'   \item BDI.       Beck Depression Inventory (BDI) score
#'   \item med.       Antidepressant medication use
#'   \item med.time.  The starting time of medication
#' }
#'
#' @docType data
#' @keywords datasets
#' @name BDIdata
#' @usage data(BDIdata)
#' @format A data frame with 7117 rows and 5 variables.
#' @references{ \enumerate{
#' \item   Wu, C. O., Tian, X. and Bang, H. A varying-coefficient model for the
#'         evaluation of time-varying concomitant intervention effects in longitudinal
#'          studies. Statistics in Medicine, 27:3042-3056, 2008.
#' \item    Wu, C. O., Tian, X. and Jiang, W. A shared parameter model for the
#'   estimation of longitudinal concomitant intervention effects Biostatistics,
#'   12(4):737-749, 2011.  }
#' }
"BDIdata"
