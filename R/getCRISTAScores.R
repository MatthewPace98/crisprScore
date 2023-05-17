#' @title ...
#' @description ...
#' 
#' @param ...
#' 
#' @details ...
#' 
#' @return A data.frame with predicted off-target score (between 0 and 1).
#'     A higher score indicates a lesser chance and severity of off-target.
#' 
#' @references ...
#' 
#' @author Matthew Pace
#' 
#' @examples ...
#' @export
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun


getCRISTAScores <- function(protospacer, spacer){
  protospacer <- .checkSequenceInputs(protospacer)
  spacer <- .checkSequenceInputs(spacer)
  if (unique(nchar(protospacer)) != 20) {
    stop("Protospacer must have length 20nt.")
  }

  if (unique(nchar(spacer)) != 29) {
    stop("Spacer must have length 29nt ([3nt][20nt][NGG][3nt]).")
  }
  results <- basiliskRun(
    env = env_crista,
    shared = FALSE,
    fork = FALSE,
    fun = .crista_python,
    protospacer = protospacer,
    spacer = spacer
  )
  return(results)
}

#' @import crisprScoreData
.crista_python <- function(protospacer, spacer){
  program <- system.file("python",
                         "crista",
                         "CRISTA.py",
                         package="crisprScore",
                         mustWork=TRUE)
  wd <- system.file("python",
                    "crista",
                    package="crisprScore")
  setwd(wd)
  df <- data.frame(sequence=protospacer,
                   score=NA_real_,
                   stringsAsFactors=FALSE)
  good <- !grepl("N", protospacer)
  protospacer.valid <- protospacer[good]
  good <- !grepl("N", spacer)
  spacer.valid <- spacer[good]
  if (length(protospacer.valid)>0){
    scores <- rep(NA_real_, length(protospacer.valid))
    for (i in seq_along(protospacer.valid)){
      protospacer <- protospacer.valid[i]
      spacer <- spacer.valid[i]
      cmd <- paste0("python ", program, " -s ", protospacer, " -d ", spacer)
      # Run the command and capture the output
      output <- system(cmd, intern = TRUE)
      score <- paste(output, collapse = "\n")
      match <- regexpr("(?<=CRISTA predicted score: )[0-9.]+", output[4], perl = TRUE)
      score <- as.numeric(regmatches(output[4], match))
      scores[i] <- score
    }
    df$score[good] <- scores
  }
  return(df)
}
