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
    stringsAsFactors = FALSE
  )
  rlt
}