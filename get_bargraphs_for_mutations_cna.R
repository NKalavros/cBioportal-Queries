bargraphs_mutations_cna <- function(gene,df){
  plotmutations <- subset(df,
                          df$Mutation.Frequency != 0)
  plotmutations$Mutation.Frequency <- plotmutations$Mutation.Frequency*100
    
  rotate_x_mut <- function(data, column_to_plot, labels_vec, rot_angle, spacebetweenbars, gene) {
    op <- par(mar=c(15,6.5,6.5,2))  
    plt <- barplot(data[[column_to_plot]], col='steelblue', xaxt="n",
                     space = rep(spacebetweenbars,length(data[[column_to_plot]])),
                     main = paste("Mutation Frequency of", gene),
                     ylab = "Percent Mutated")
      text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1,1), xpd = TRUE, cex=0.6) 
  }
  pdf(paste(gene,"MutFreq.pdf"), width = 15, height = 10)
  rotate_x(plotmutations,"Mutation.Frequency", rownames(plotmutations),45, 0.25, gene)
  dev.off()
  
  plotcna <- subset(df,
                    (df$Gain.Frequency != 0 | !is.na(df$Gain.Frequency)) | (df$Amplification.Frequency != 0 | !is.na(df$Amplification.Frequency)) | (df$Shallow.loss.Frequency != 0 | !is.na(df$Shallow.loss.Frequency)) | (df$Deep.deletion.Frequency != 0 | !is.na(df$Deep.deletion.Frequency)))
  
  plotcna <- plotcna*100
  
  rotate_x_cna <- function(data, column_to_plot, labels_vec, rot_angle, spacebetweenbars, gene) {
    op <- par(mfrow=c(1, 1), mar=c(15,7.5,7.5,7.5))  
    plt <- barplot(t(as.matrix(data[column_to_plot])), col= c('steelblue',"#727272","#f1595f", "#79c36a"), xaxt="n",
                   space = rep(spacebetweenbars,ncol(data[column_to_plot])),
                   main = paste("Copy Number Alteration Frequency of", gene),
                   ylab = "Percent Mutated",
                   legend = c("Gain Frequency","Amplification Frequency","Shallow Loss Frequency","Deep Deletion Frequency"),
                   args.legend = list(x = ncol(t(as.matrix(data[column_to_plot]))) + 26.5,
                                      y = max(colSums(t(as.matrix(data[column_to_plot])))) + 10,
                                      bty = "n"))
    text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1,1), xpd = TRUE, cex=0.6) 
  }
  pdf(paste(gene,"CnaFreq.pdf"), width = 15, height = 10)
  
  rotate_x_cna(plotcna,colnames(plotcna)[2:5], rownames(plotcna),45, 0.25, gene)
  dev.off()
}