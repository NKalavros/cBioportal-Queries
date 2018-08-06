#Function accepts a gene name in Onco Query language and a boolean to check whether or not to query provisional studies
obtain_study_data_mutations_cna <- function(genes,queryprovisionals){
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
  geneticprofilesidsmutations <- rep(NA,length(data))
  
  #For every data, if there is a category in its genetic_profile_id
  #Named mutations, then save that specific name in the empty vector
  #Created earlier
  
  for(i in seq(1:length(data))) {
    if (sum(grepl("mutations",data[[i]][["genetic_profile_id"]],
                  ignore.case = TRUE)) > 0) {
      geneticprofilesidsmutations[i] <- data[[i]][["genetic_profile_id"]][grepl("mutation",
                                                                                data[[i]][["genetic_profile_id"]],
                                                                                ignore.case = TRUE)]
    }
  }              
  
  #Same for gistic
  geneticprofileidscna <- rep(NA,length(data))
  
  for(i in seq(1:length(data))) {
    if (sum(grepl("gistic",data[[i]][["genetic_profile_id"]],
                  ignore.case = TRUE)) > 0) {
      geneticprofileidscna[i] <- data[[i]][["genetic_profile_id"]][grepl("gistic",
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
  caselistsmutations <- rep(NA,length(data))
  
  caselistscna <- rep(NA, length(data))
  
  #Get the caselist name that has all the sequenced tumors.
  for(i in seq(1:length(datacaselists))){
    if(sum(grepl("sequenced",datacaselists[[i]],
                 ignore.case = TRUE)) > 0){
      caselistsmutations[i] <- datacaselists[[i]][grepl("sequenced",datacaselists[[i]],
                                                      ignore.case = TRUE)]
    }
  }
  
  for(i in seq(1:length(datacaselists))){
    if(sum(grepl("cna$",datacaselists[[i]],
                 ignore.case = TRUE)) > 0){
      caselistscna[i] <- (datacaselists[[i]][grepl("cna$",datacaselists[[i]],
                                                   ignore.case = TRUE)])
    
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
  
  mygenemut <- mapply(foo3,geneticprofilesidsmutations,
                   caselistsmutations,
                   SIMPLIFY = TRUE)
  
  #Create an empty vector
  mutfreq <- rep(NA,length(mygenemut))
  
  #Find frequencies
  for(i in seq(1:length(mygenemut))){
    mutfreq[i] <- sum(!is.na(mugenemut[[i]]))/length(mygenemut[[i]])
    
  }
  names(mutfreq) <- queriedstudies[,2]
  
  #Same thing for the CNAs
  mygenecna <- mapply(foo3,geneticprofileidscna,caselistscna,
                   SIMPLIFY = TRUE)
  
  #Create 4 empty vectors to represent shallow deletion, deep deletion, gain and amplification.
  gainfreq <- rep(NA, length(mygenecna))
  doublegainfreq <- rep(NA,length(mygenecna))
  shallowlossfreq <- rep(NA,length(mygenecna))
  deepdeletionfreq <- rep(NA,length(mygenecna))
  
  for(i in seq(1:length(mygenecna))){
    gainfreq[i] <- sum(mygenecna[[i]] == 1)/length(mygenecna[[i]])
    doublegainfreq[i] <- sum(mygenecna[[i]] == 2)/length(mygenecna[[i]])
    shallowlossfreq[i] <- sum(mygenecna[[i]] == -1)/length(mygenecna[[i]])
    deepdeletionfreq[i] <- sum(mygenecna[[i]] == -2)/length(mygenecna[[i]])
  }
  
  names(gainfreq) <- names(doublegainfreq) <- names(shallowlossfreq) <- names(deepdeletionfreq) <- queriedstudies[,2]
  
  result <- list(mutfreq,gainfreq,doublegainfreq,shallowlossfreq,deepdeletionfreq)
  
  names(result) <- c("Mutation Frequency","Gain Frequency","Amplification Frequency","Shallow loss Frequency","Deep deletion Frequency")
  
  return(result)
}




pnpla2 <- obtain_study_data_mutations_cna("PNPLA2", TRUE)
lipe <- obtain_study_data_mutations_cna("LIPE", TRUE)
