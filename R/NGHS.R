#' NGHS dataset
#'
#' This dataset includes 2378  girls (total 19701 observations) enrolled in the
#' National Heart, Lung, and Blood Institute's Growth and Health Study (NGHS).
#' NGHS is a multicenter population-based cohort study aimed at evaluating
#' the racial differences and longitudinal changes in childhood cardiovascular
#' risk factors between  Caucasian and  African American girls during
#' childhood and adolescence.
#'
#' \itemize{
#'   \item ID.        Subject ID
#'   \item RACE.     Subject's race (1=Caucasian, 2= African American)
#'   \item AGE,HEIGHT,WEIGHT,BMI.       Age, height, weight and BMI at study visit
#'   \item BMIPCT, HTPCT.    CDC Age-adjusted BMI percentile and height percentile at study visit
#'   \item SBP,DBP.    Systolic and diastolic blood pressure at study visit
#'   \item TG,LDL.     Triglyceride and Low-density lipoprotein (LDL) cholesterol at study visit
#' }
#'
#' @docType data
#' @keywords datasets
#' @name NGHS
#' @usage data(NGHS)
#' @format A data frame with 19701 rows and 12 variables.
#' @references{ \enumerate{
#' \item   National Heart, Lung, and Blood Institute Growth and Health Research Group (NGHSRG).
#'          Obesity and cardiovascular disease risk factors in black and white girls: the NHLBI Growth and Health Study.
#'          American Journal of Public Health, 82:1613-1620, 1992.
#' \item    Wu, C. O. and Tian, X. Nonparametric estimation of conditional distributions
#' and rank-tracking probabilities with time-varying transformation
#' models in longitudinal studies. Journal of the American Statistical Association,
#' 108:971-982, 2013.}
#' }
"NGHS"


