## clean data and extract variables we want to analyze
## use the functions in cleandata.R
# Wed Oct 31 18:12:27 2018 ------------------------------


# function of getting cognitive function scores
dir <- "~/AsynchronousLongitudinal/RealData/"

# Cognitive functions
# folders that have cognitive function scores
folders <- c("ICPSR_V4", "ICPSR_V6", "ICPSR_V7", "ICPSR_V8", "ICPSR_V9", "ICPSR_V10")
cog.vars <- getVariables(dir, folders, vars = c("IMED", "SDMT", "DIGIT", "DLAY"), type = "COGDAY", strict = FALSE)
# remove observations with no cognitive function days
#cog.vars <- cog.vars[!is.na(cog.vars$COGDAY), ]
cog.data <- getCogScore(cog.vars)

# Physical Measures
folders <- list.files("../")[grep("ICPSR", list.files("../"))]
phy.data <- getVariables(dir, folders, vars = c("PULSE","SYSBP1","DIABP1","SYSBP2","DIABP2","BMI"),
                         type = "PHYDAY", strict = TRUE)
phy.data <- phy.data %>% 
           mutate(DIABP = (DIABP1+DIABP2)/2,
                  SYSBP = (SYSBP1+SYSBP2)/2) %>% as.data.frame
for (var in c("DIABP1", "DIABP2", "SYSBP1", "SYSBP2")) {phy.data[[var]] <- NULL}

# Harmone Measures, don not include TSH since samples are not enough
hrm.data <- getVariables(dir, folders, vars = c("SHBG","DHAS", "FSH"), type = "HRMDAY", strict = TRUE)

# Cardiovascular Measures
cvr.data <- getVariables(dir, folders, vars = c("CHOLRES","GLUCRES","INSURES", "TRIGRES"), type = "CVRDAY", strict = TRUE)

# join the three measures, 11 X
full.data <- full_join(cog.data, phy.data, by = c("ID", "VISITS")) %>%
            full_join(., hrm.data, by = c("ID", "VISITS")) %>%
            full_join(., cvr.data, by = c("ID", "VISITS")) %>% as.data.frame
full.data <- full.data[order(full.data$ID, full.data$VISITS), ]
rownames(full.data) <- NULL
saveRDS(full.data, file = "./CognitiveScores_Xvars.rds")
write.csv(full.data, file = "./CognitiveScores_Xvars.csv")

#########################################################
#####################Z variables#########################
#########################################################
# variables from baseline dataset
folder <- "ICPSR_V0"
files <- list.files(paste0(dir, folder, "/DS0001"))
file <- files[grep("rda", files)]
data <- get(load(paste0(dir, folder, "/DS0001/", file)))
## variables extracting
vars.str <- paste(c("SWANID", "\\bAGE0\\b", "OCCUP0", "RACE", "RELIPRE0", "HOUSEHL0"), collapse = "|")
z.data <- z_process(vars.str, data, decode = FALSE)
z.data$OCCUP0 <- unlist(sapply(z.data$OCCUP0, rsps_coding))
z.data$RACE <- unlist(sapply(z.data$RACE, rsps_coding))
z.data$RELIPRE0 <- unlist(sapply(z.data$RELIPRE0, rsps_coding))
z.data$HOUSEHL0 <- unlist(sapply(z.data$HOUSEHL0, rsps_coding))
### smoking variables
z.data <- cbind(z.data, z_process("STRTSMO|AVCIGDA", data, decode = FALSE))
### working variables
# work hours per day, EVENING/SWING (Between 3 PM and 11 PM) and NIGHT (Between 9 PM and 9 AM)
z.data <- cbind(z.data, z_process("DAYSHFT|EVESHFT|NGHTSHF|ROTSHFT", data))
## working movement
z.data <- cbind(z.data, z_process("SIT0|STAND0|WALK0|LIFT0|STOOP0|PUSH0|SWEAT0", data))
### concerning language
z.data$LANGREA <- unlist(sapply(data$LANGREA0, rsps_coding)) # language you speak/read
z.data$LANGTHN <- unlist(sapply(data$LANGTHN0, rsps_coding)) # language you think
# bilingual indicator
z.data$BILINGUAL <- sapply(as.numeric(z.data$LANGREA)-as.numeric(z.data$LANGTHN), function(x) {ifelse(x!=0, 1, 0)})
### personal feelings
z.data <- cbind(z.data, 
                z_process("COURTES0|RESPECT0|POORSER0|NOTSMAR0|AFRAIDO0|DISHONS0|BETTER0|INSULTE0|HARASSE0|IGNORED0", data))
### personal assistant
z.data <- cbind(z.data, z_process("LISTEN0|TAKETOM0|CONFIDE0|HELPSIC0", data))
## how you pay the insurance
z.data <- cbind(z.data, z_process("PREPAID0|OTHRPRI0|MEDICAR0|MEDICAI0|MILITAR0|NOINSUR0|OTHINSU0", data))

## remove the "0" at the end of variable names
names(z.data) <- sapply(names(z.data), function(x) {
  if (substr(x, nchar(x), nchar(x)) == "0") {split_str(x, last = FALSE)}
  else {x}
})
saveRDS(z.data, file = "./Zvars.rds")
write.csv(z.data, file = "./Zvars.csv")


###################################################################
###################Nutrition Variables#############################
###################################################################
folders <- c("ICPSR_V0", "ICPSR_V5", "ICPSR_V9")
dir <- "~/AsynchronousLongitudinal/RealData/"
nutr.vars <- c("ALL", "DTTNIAC")
#vars = c("DTT", "GMSOLID", "FIBBEAN", "FIBVEGF", "FIBGRAI")
nutr.data1 <- getVariables(dir, folders[1], vars = nutr.vars,
                         type = "FFQDAY", strict = FALSE)
nutr.data5 <- getVariables(dir, folders[2], vars = nutr.vars,
                           type = "FFQDAY", strict = FALSE)
nutr.data9 <- getVariables(dir, folders[3], vars = nutr.vars,
                           type = "FFQDAY", strict = FALSE)
com.vars <- intersect(intersect(names(nutr.data1), names(nutr.data5)), names(nutr.data9))
nutr.data <- rbind(nutr.data1[, com.vars], nutr.data5[, com.vars], nutr.data9[, com.vars])
nutr.data <- nutr.data[order(nutr.data$ID, nutr.data$VISITS), ]
rownames(nutr.data) <- NULL

saveRDS(nutr.data, file = "./NUTRITIONvars.rds")
write.csv(nutr.data, file = "./NUTRITIONvars.csv")
