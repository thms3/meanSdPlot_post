#!/usr/bin/env Rscript
#@uthor : Thomas NEFF

# lib & src ---------------------------------------------------------------
library(MASS)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ggnewscale)

# data --------------------------------------------------------------------
load(file = "/home/thomas/Documents/PNC_PNET/rdata/pair_idhmut/pair_idhmut_dds.RData")

# exe ---------------------------------------------------------------------
counts <- DESeq2::counts(dds_wald$data)
counts <- DESeq2::rlog(object = counts,blind = F)

res <- vsn::meanSdPlot(counts,plot = F,ranks = F)

meanSDplot <- function(counts,genes=NULL,results=NULL,pval=0.05,colname_pval="padj",colname_fold="log2FoldChange"
                       ,rank=F,n_density=100,sel=list(size=2,shape=21,colors=c("red","grey"),labels=c("up","down"),name="Exp."),
                       show.labels=T,graph=T,alpha_dens=1){

  mean <- apply(X = counts,MARGIN = 1,FUN = function(x) mean(x))
  sd <- apply(X = counts, MARGIN = 1, FUN = function(x) sd(x))
  df <- data.frame(merge(x = data.frame(mean),y = data.frame(sd),by="row.names",full=T),row.names = "Row.names")
  
  density <- MASS::kde2d(x = df$mean, y = df$sd,n=n_density)
  xi <- findInterval(df$mean, density$x)
  yi <- findInterval(df$sd, density$y)
  df$density <- density$z[cbind(xi,yi)] * length(df$mean)/n_density
  
  df <- df[order(df$mean),]
  df$rank <- seq(1,length(df$mean))
  
  dm = 0.025
  quant <- seq(dm, 1-dm, by = dm)
  inframe <- function(x, x1, x2) { (x >= x1) & (x <= x2) }
  rank.sd <- sapply(quant, function(y) median(df$sd[inframe(x = df$rank/length(df$rank), x1 = y - 2*dm, x2 = y + 2*dm)], na.rm = TRUE)) 
  rank.quant <- seq(1,length(df$mean),length.out = 41)[-c(1,41)]
  quantile <- sapply(rank.quant, function(x) df$mean[round(x)])
  
  if (!rank) { 
    df$rank <- df$mean
    rank.quant <- quantile
    }

  line.sd.quant <- data.frame(cbind(quantile,rank.sd,rank.quant))
  
  if (!is.null(results)) {
    results$test <- ifelse(results[,colname_pval] > pval | is.na(results[,colname_pval]),"no.sign","sign")
    results$dir <- ifelse(sign(results[,colname_fold])==1,"up","down")
    df <- data.frame(merge(x=df,y = data.frame(results[,c("test","dir")]),by="row.names",full=T),row.names="Row.names")
  }
  
  if (!is.null(genes)) {
    df$test <- sapply(rowmames(df),function(x) ifelse(test = x %in% genes,"sign","no.sign"))
  }
  
  if (is.null(genes) & is.null(results)) {
    p <-  ggplot() + geom_point(data = df, mapping = aes(x = rank, y = sd, color = density),alpha=alpha_dens) +
      scale_color_continuous(name="counts") + ggnewscale::new_scale_color() +
      geom_line(data = line.sd.quant, mapping = aes(x = rank.quant,y = rank.sd,color="median")) + 
      scale_color_manual(name="",values=c("red"))
  } else {
    p <- ggplot() + geom_point(data = df[df$test=="no.sign",], mapping = aes(x = rank, y = sd, color = density)) +
      geom_point(data = df[df$test=="sign",],mapping = aes(x = rank,y = sd,fill=dir),size=sel$size,shape=sel$shape) + 
      scale_fill_manual(name=sel$name,values = sel$colors,labels=sel$labels,aesthetics = "fill") +
      scale_color_continuous(name="counts") + ggnewscale::new_scale_color() +
      geom_line(data = line.sd.quant, mapping = aes(x = rank.quant,y = rank.sd,color="median")) + 
      scale_color_manual(name="",values=c("red"))
  }
  
  if (show.labels & (!is.null(genes) | !is.null(results))) {
    p <- p + ggrepel::geom_text_repel(data = df[df$test=="sign",],mapping = aes(x = rank,y = sd),label=rownames(df[df$test=="sign",]))
  }
  if (graph) {  print(p) }


  return(df)
}

df_genes <- meanSDplot(counts,results = dds_wald$res,rank = F)


  