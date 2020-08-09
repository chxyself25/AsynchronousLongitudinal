## clean data functions
## we need 10 X(time-variant) variables, 10 Z(time-invariant) variables
# Wed Oct 31 17:02:26 2018 ------------------------------
library(dplyr)
library(ggplot2)
library(readr)

# function for detecting response code
# response 1 as 1, the others are zero
rsps_code <- function(y) {
  if (is.na(y)) {
    return(NA)
  }
  else {
    if (grepl("1", y)) {
      return(1)
    }
    else {
      return(0)
    }
  }
}
# function that remove the brackets and return the number inside
rsps_coding <- function(y) {
  if(is.na(y)) {
    return(NA)
  } else {
    y.str <- strsplit(as.character(y), split = " ")[[1]]
    str.len <- length(strsplit(y.str[1], split = "")[[1]])
    return(substr(y.str[1], 2, (str.len-1)))
  }
}
## function that is able to remove the number at the end of variable names
split_str <- function(x, last) {
  split.str <- strsplit(x, split = "")[[1]]
  n <- length(split.str)
  if (last) {
    return(paste(split.str[1:(n-2)], collapse = ""))
  } else {
    return(paste(split.str[1:(n-1)], collapse = ""))
  }
}

# function that extracts desired variables
# vars is variable characters, type is the observation days
# folders are diectories, dir is directory where these folders stored(should have /)
getVariables <- function(dir, folders, vars, type, strict = FALSE) {
  out.data <- NULL
  for (folder in folders) {
    files <- list.files(paste0(dir, folder, "/DS0001/"))
    file <- files[grep(".tsv", files)]
    if (length(file) != 0) {
      data <- read_delim(paste0(dir, folder, "/DS0001/", file), "\t", escape_double = FALSE, trim_ws = TRUE) %>% as.data.frame
    } else {
      file <- files[grep("rda", files)]
      data <- get(load(paste0(dir, folder, "/DS0001/", file)))
    }
    visit.no <- gsub("ICPSR_V", "", folder)
    if (strict) {
      var.str <- paste(paste("\\b", c(type, vars), visit.no, "\\b", sep = ""), collapse = "|")
    } else {
      var.str <- paste(c(type, vars), collapse = "|")
    }
    day.idx <- names(data)[grep(var.str, names(data))] %>% sort
    subdata <- data[, day.idx]
    # if there is no observation time, skip this folder
    if (all(!grepl(type, names(subdata)))) {
      next
    }
    # remove the numbers in variable names
    names(subdata) <- sapply(names(subdata), split_str, last = grepl("10", folder))
    subdata$ID <- data[,1] # the first column is ID
    subdata$VISITS <- as.numeric(visit.no)
    if ("SDMT" %in% vars) {subdata$SDMTSPE <- NULL} # some has SDMTSPE some not
    if ("BMI" %in% vars) {
      if (all(!grepl("BMI", names(subdata)))) {
        subdata$BMI <- NA
      }
    }
    # if ("TSH" %in% vars) {
    #   if (all(!grepl("TSH", names(subdata)))) {
    #     subdata$TSH <- NA
    #   }
    # }
    # remove observations with no days
    #out.data <- rbind(out.data, subdata[!is.na(subdata[[type]]), ])
    out.data <- rbind(out.data, subdata)
  }
  out.data$ID <- as.character(out.data$ID)
  return(out.data)
}

## function that calculate cognitive function scores by using different methods
getCogScore <- function(whdata) {
  whdata <- whdata[order(whdata$ID, whdata$COGDAY), ]
  var_names <- names(whdata)
  # DIGIT BACKWARDS SCORE
  digits_back <- whdata[, grep("DIGIT", var_names)]
  digits_back <- apply(digits_back, 2, function(x) {unlist(sapply(x, rsps_code))})
  whdata$DIGIT <- apply(digits_back, 1, sum, na.rm = TRUE)
  # DELAYED RECALL OF STORY SCORE
  delay <- whdata[, grep("DLAY", var_names)]
  delay <- apply(delay, 2, function(x) {unlist(sapply(x, rsps_code))})
  whdata$DELAY <- apply(delay, 1, sum, na.rm = TRUE)
  # IMMEDIATE RECALL OF STORY SCORE
  immediate <- whdata[, grep("IMED", var_names)]
  immediate <- apply(immediate, 2, function(x) {unlist(sapply(x, rsps_code))})
  whdata$IMED <- apply(immediate, 1, sum, na.rm = TRUE)
  # SDMT score
  whdata$SDMT <- (whdata$SDMTCOR/whdata$SDMTATM)*12 # to make them have the same scale
  # z scores
  whdata$digit.z <- with(whdata, (DIGIT-mean(DIGIT, na.rm = TRUE))/sd(DIGIT, na.rm = TRUE))
  whdata$delay.z <- with(whdata, (DELAY-mean(DELAY, na.rm = TRUE))/sd(DELAY, na.rm = TRUE))
  whdata$imed.z <- with(whdata, (IMED-mean(IMED, na.rm = TRUE))/sd(IMED, na.rm = TRUE))
  whdata$sdmt.z <- with(whdata, (SDMT-mean(SDMT, na.rm = TRUE))/sd(SDMT, na.rm = TRUE))
  # final score for cognitive function performance
  whdata <- whdata %>% 
    mutate(COGSCORE = DIGIT + DELAY + IMED + SDMT,
           GLBSCORE = (digit.z + delay.z + imed.z + sdmt.z)/4)
  cog.data <- whdata[, c("ID", "COGDAY", "COGSCORE", "GLBSCORE", "IMED", "DELAY", "DIGIT", 
                         "SDMT", "SDMTCOR", "SDMTATM", "VISITS")]
  return(cog.data)
}

## function that processes z variables: time-invariant variables
z_process <- function(z_str, data, decode = TRUE){
  var.idx <- names(data)[grep(z_str, names(data))] %>% sort
  z.data <- data[, var.idx]
  if (decode) {
    for (var in var.idx) {
      z.data[[var]] <- unlist(sapply(z.data[[var]], rsps_coding))
    } 
  }
  return(z.data)
}

