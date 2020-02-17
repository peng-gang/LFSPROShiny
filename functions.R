load("data/mortality_incidence_all_cancer.RData")

format.pv <- function(pv){
  if(pv < 0.001){
    return("<0.001")
  } else {
    return(format(round(pv, digits = 3), nsmall=3))
  }
}

lfs.match <- function(n.total, lrt, idx){
  rlt <- rep(NA, n.total)
  rlt[idx] <- lrt
  rlt
}


runLFSPRO <- function(fam.data, cancer.data, counselee.id){
  rlt.chompret <- lfsChompret2015(fam.data, cancer.data, counselee.id)
  rlt.classic <- lfsClassic(fam.data, cancer.data, counselee.id)
  rlt.lfspro <- lfspro(fam.data, cancer.data, counselee.id)
  
  nc.5 <- NULL
  nc.10 <- NULL
  nc.15 <- NULL
  for(i in 1:nrow(counselee.id)){
    idx <- which(fam.data$id == counselee.id$id[i])
    #print(idx)
    #print(fam.data$age[idx])
    #print(ifelse(fam.data$gender[idx]==0, "female", "male"))
    nc.5 = c(nc.5, riskCauseSpecific(fam.data$age[idx], fam.data$age[idx] + 5, ifelse(fam.data$gender[idx]==0, "female", "male"), incidence_all_cancer, mortality))
    nc.10 = c(nc.10, riskCauseSpecific(fam.data$age[idx], fam.data$age[idx] + 10, ifelse(fam.data$gender[idx]==0, "female", "male"), incidence_all_cancer, mortality))
    nc.15 = c(nc.15, riskCauseSpecific(fam.data$age[idx], fam.data$age[idx] + 15, ifelse(fam.data$gender[idx]==0, "female", "male"), incidence_all_cancer, mortality))
  }
  
  idx.cs <- match(rlt.lfspro$Cancer_specific_risks$Breast_risks$counselee.id, counselee.id$id)
  idx.mpc <- match(rlt.lfspro$Multiple_primary_cancer_risks$id, counselee.id$id)
  
  rlt <- data.frame(
    id = counselee.id$id,
    classic = rlt.classic$result,
    chompret = rlt.chompret$result,
    carrier = rlt.lfspro$Mutation_probability$mutation_probability,
    breast.5 = lfs.match(length(counselee.id$id), 
                        as.numeric(as.character(rlt.lfspro$Cancer_specific_risks$Breast_risks$`breast in 5 yrs`)), 
                        idx.cs),
    breast.10 = lfs.match(length(counselee.id$id), 
                          as.numeric(as.character(rlt.lfspro$Cancer_specific_risks$Breast_risks$`breast in 10 yrs`)), 
                          idx.cs),
    breast.15 = lfs.match(length(counselee.id$id), 
                          as.numeric(as.character(rlt.lfspro$Cancer_specific_risks$Breast_risks$`breast in 15 yrs`)), 
                          idx.cs),
    sarcoma.5 = lfs.match(length(counselee.id$id),
                          as.numeric(as.character(rlt.lfspro$Cancer_specific_risks$Sarcoma_risks$`sarcoma in 5 yrs`)),
                          idx.cs),
    sarcoma.10 = lfs.match(length(counselee.id$id),
                          as.numeric(as.character(rlt.lfspro$Cancer_specific_risks$Sarcoma_risks$`sarcoma in 10 yrs`)),
                          idx.cs),
    sarcoma.15 = lfs.match(length(counselee.id$id),
                          as.numeric(as.character(rlt.lfspro$Cancer_specific_risks$Sarcoma_risks$`sarcoma in 15 yrs`)),
                          idx.cs),
    other.5 = lfs.match(length(counselee.id$id),
                        as.numeric(as.character(rlt.lfspro$Cancer_specific_risks$Other_cancers$`others in 5 yrs`)),
                        idx.cs),
    other.10 = lfs.match(length(counselee.id$id),
                        as.numeric(as.character(rlt.lfspro$Cancer_specific_risks$Other_cancers$`others in 10 yrs`)),
                        idx.cs),
    other.15 = lfs.match(length(counselee.id$id),
                        as.numeric(as.character(rlt.lfspro$Cancer_specific_risks$Other_cancers$`others in 15 yrs`)),
                        idx.cs),
    second.5 = lfs.match(length(counselee.id$id),
                         rlt.lfspro$Multiple_primary_cancer_risks$`5 years`, idx.mpc),
    second.10 = lfs.match(length(counselee.id$id),
                         rlt.lfspro$Multiple_primary_cancer_risks$`10 years`, idx.mpc),
    second.15 = lfs.match(length(counselee.id$id),
                         rlt.lfspro$Multiple_primary_cancer_risks$`15 years`, idx.mpc),
    nc.5 = nc.5,
    nc.10 = nc.10,
    nc.15 = nc.15,
    stringsAsFactors = FALSE
  )
  rlt
}



riskNoCompeting <- function(age1, age2, sex = c("male", "female"), incidence_all_cancer){
  # cancer risk between age1 to age2
  # age1 <= T <= age2
  sex <- match.arg(sex)
  if(sex == "male"){
    S1 <- exp(-sum(incidence_all_cancer$male[1:age1])) 
    S2 <- exp(-sum(incidence_all_cancer$male[1:(age2+1)]))
    rlt <- (S1 - S2)/S1
  } else if(sex == "female") {
    S1 <- exp(-sum(incidence_all_cancer$female[1:age1])) 
    S2 <- exp(-sum(incidence_all_cancer$female[1:(age2+1)]))
    rlt <- (S1 - S2)/S1
  } else {
    stop("Gender should be male or female")
  }
  rlt
}


riskCauseSpecific <- function(age1, age2, sex = c("male", "female"), incidence_all_cancer, mortality){
  # cancer risk between age1 to age2
  # age1 <= T <= age2
  sex <- match.arg(sex)
  rlt <- NULL
  if(sex == "male"){
    S1 <- exp(-sum(incidence_all_cancer$male[1:age1])) 
    riskAge <- 0
    for(i in age1:age2){
      riskAge <- riskAge + riskCauseSpecificAge(i, "male", incidence_all_cancer, mortality)
    }
    rlt <- riskAge/S1
  } else if(sex == "female"){
    S1 <- exp(-sum(incidence_all_cancer$female[1:age1])) 
    riskAge <- 0
    for(i in age1:age2){
      riskAge <- riskAge + riskCauseSpecificAge(i, "female", incidence_all_cancer, mortality)
    }
    rlt <- riskAge/S1
  } else {
    stop("Gender should be male or female")
  }
  rlt
}


riskCauseSpecificAge <- function(age, sex = c("male", "female"), incidence_all_cancer, mortality){
  # cancer risk at age
  rlt <- 0
  sex <- match.arg(sex)
  
  if(sex == "male"){
    if(age==0){
      rlt <- incidence_all_cancer$male[1]
    } else if(age >=85 ) {
      stop("Age should be less than 85")
    } else {
      S = exp(-(sum(incidence_all_cancer$male[1:age]) + sum(mortality$male[1:age])))
      rlt <- S * incidence_all_cancer$male[age + 1]
    }
  } else if(sex == "female"){
    if(age==0){
      rlt <- incidence_all_cancer$female[1]
    } else if(age >=85 ) {
      stop("Age should be less than 85")
    } else {
      S = exp(-(sum(incidence_all_cancer$female[1:age]) + sum(mortality$female[1:age])))
      rlt <- S * incidence_all_cancer$female[age + 1]
    }
  } else {
    stop("Gender should be male or female")
  }
  rlt
}


