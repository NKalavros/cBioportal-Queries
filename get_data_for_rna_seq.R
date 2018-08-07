#Takes the same input as others plus a z-score threshold (positive integer)
#Returns frequency of alterations in all queried studies.
obtain_study_data_rnaseq <- function(genes,queryprovisionals,threshold){
  #Load the package
  require(cgdsr)
  
  #Create a connection object
  
  mycgds <- CGDS("http://www.cbioportal.org/")
  
  #Test it
  
  test(mycgds)
  
  #Get all studies
  studies <- getCancerStudies(mycgds)
  
  #Which aren't pancancer?
  nopancancer <- !grepl("pan_can",studies[,c(1)])
  
  #Which aren't provisional? (For now let's keep the provisionals)
  if(queryprovisionals != TRUE){
    noprovisional <- !grepl("Provisional", studies[,c(2)])
  }else{
    noprovisional <- rep(TRUE,length(studies[,c(1)]))
  }
  
  #Keep only them, or not the provisionals, depending on the boolean value you entered
  newstudies <- nopancancer & noprovisional
  
  #Get the studies' names
  queriedstudies <- studies[newstudies,]
  
  #Place holder function because I hate passing function arguments to 
  #lapply
  foo <- function(studyname){
    return(getGeneticProfiles(mycgds,studyname))
  }
  
  #Get the genetic profiles for each of those studies
  data <- lapply(queriedstudies[,1], FUN = foo)
  
  #Create an empty vector of length the same as your data
  geneticprofilesidszscores <- rep(NA,length(data))
  
  #For every data, if there is a category in its genetic_profile_id
  #Named mutations, then save that specific name in the empty vector
  #Created earlier
  
  for(i in seq(1:length(data))) {
    if (sum(grepl("zscores",data[[i]][["genetic_profile_id"]],
                  ignore.case = TRUE)) > 0) {
      geneticprofilesidszscores[i] <- data[[i]][["genetic_profile_id"]][grepl("zscores",
                                                                                data[[i]][["genetic_profile_id"]],
                                                                                ignore.case = TRUE)]
    }
  }              
  
  #Same placeholder function as before
  foo2 <- function(studyname){
    return(getCaseLists(mycgds,studyname)[,c(1)])
  }
  
  #Get all the caselists
  datacaselists <- lapply(queriedstudies[,c(1)], foo2)
  
  #Create empty vector
  caselistzscores <- rep(NA,length(data))
  
  #Get the caselist name that has all the sequenced tumors.
  for(i in seq(1:length(datacaselists))){
    if(sum(grepl("rna_seq",datacaselists[[i]],
                 ignore.case = TRUE)) > 0){
      caselistzscores[i] <- datacaselists[[i]][grepl("rna_seq",datacaselists[[i]],
                                                        ignore.case = TRUE)]
    }
  }
  
  #Same thing with the placeholder functions, I will make them
  #Generic after the 4th I guess, it is just not worth it
  foo3 <- function(geneticprofilesid,caselist){
    if((is.na(geneticprofilesid) | is.na(caselist)) == TRUE){
      return(NA)
    }else{
      return(c(getProfileData(x = mycgds,
                              genes = genes,
                              geneticProfiles = geneticprofilesid,
                              caseList = caselist)))
    }
  }
  
  pnpla2rna <- mapply(foo3,geneticprofilesidszscores,
                      caselistzscores,
                      SIMPLIFY = TRUE)
  
  #Create 4 empty vectors to represent shallow deletion, deep deletion, gain and amplification.
  upregfreq <- rep(NA, length(pnpla2rna))
  downregfreq <- rep(NA,length(pnpla2rna))
  
  for(i in seq(1:length(pnpla2rna))){
    upregfreq[i] <- sum(pnpla2rna[[i]] > threshold)/length(pnpla2rna[[i]])
    downregfreq[i] <- sum(pnpla2rna[[i]] < -threshold)/length(pnpla2rna[[i]])
  }
  
  names(upregfreq) <- names(downregfreq) <- queriedstudies[,2]
  
  result <- data.frame(cbind(upregfreq,downregfreq))
  
  rownames(result) <- names(upregfreq)
  
  colnames(result) <- c("Upregulated","Downregulated")
  
  return(result)
}
