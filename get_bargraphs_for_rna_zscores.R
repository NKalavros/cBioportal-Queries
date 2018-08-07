bargraphs_rna_zscores <- function(gene,df){

  plotzscores <- subset(df,
                        df$Upregulated != 0 | df$Downregulated != 0)
  plotzscores <- plotzscores*100
  
  df <- plotzscores
  
  rotate_x_rna <- function(data, column_to_plot, labels_vec, rot_angle, spacebetweenbars, gene) {
    op <- par(mfrow=c(1, 1), mar=c(15,7.5,7.5,7.5))  
    plt <- barplot(t(as.matrix(data[column_to_plot])), col= c('#f1595f',"steelblue"), xaxt="n",
                   space = rep(spacebetweenbars,ncol(data[column_to_plot])),
                   main = paste("Percent of patients with dysregulation", gene),
                   ylab = "Percent dysregulated",
                   legend = c("Upregulated","Downregulated"),
                   args.legend = list(x = ncol(t(as.matrix(data[column_to_plot]))) + 3,
                                      y = max(colSums(t(as.matrix(data[column_to_plot])))),
                                      bty = "n"))
    text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1,1), xpd = TRUE, cex=0.6) 
  }
  pdf(paste(gene,"RnaFreq.pdf"), width = 15, height = 10)
  
  rotate_x_rna(df,colnames(df)[1:2], rownames(df),45, 0.25, gene)
  dev.off()
}
