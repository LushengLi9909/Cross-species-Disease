require(visNetwork)
library(ggplot2)
signif_results_h <- read.delim("data/signif_results_h.txt",header = TRUE,sep=",")
signif_results_m <- read.delim("data/signif_results_Mondo.txt",header = TRUE,sep=",")
signif_results_m = signif_results_m[(signif_results_m$MONDO_ID %in% signif_results_h$MONDO),]

mapping_all <- read.delim("data/ph_mp_mapping.txt",header = TRUE,sep=",")
cell_order <- read.delim("data/cell_order.txt",header = TRUE,sep=",")
net_cross_do <- function(disease){
  if(disease != ""){
    if(disease %in% signif_results_h$Disease){
      do_results_m = signif_results_m[signif_results_m$MONDO_label == disease,]
      do_results_h = signif_results_h[signif_results_h$Disease == disease,]
}else{
  disease = unique(signif_results_m[signif_results_m$DO.Disease == disease,"MONDO_label"])
  do_results_m = signif_results_m[signif_results_m$MONDO_label == disease,]
  do_results_h = signif_results_h[signif_results_h$Disease == disease,]
}

adjacency = data.frame()
adjacency = rbind(adjacency,data.frame("from"=do_results_m$Phenotype,"to"=do_results_m$CellType))
adjacency = rbind(adjacency,data.frame("from"=do_results_h$Phenotype,"to"=do_results_h$CellType))

for( i in unique(do_results_m$MPID)){
  p_h = mapping_all[mapping_all$MPID == i,"HPID"]
  if(sum((p_h %in% do_results_h$HPID) == TRUE)>0){
    adjacency = rbind(adjacency,data.frame("from"=unique(do_results_m[do_results_m$MPID ==i,"Phenotype"]),
                                           "to"=unique(do_results_h[do_results_h$HPID ==p_h[which((p_h %in% do_results_h$HPID)==TRUE)],"Phenotype"])))
  }
}
adjacency = adjacency[!duplicated(adjacency),]
pheno_m = unique(do_results_m$Phenotype)
pheno_h = unique(do_results_h$Phenotype)
Npheno_m = length(pheno_m)
Npheno_h = length(pheno_h)
cell = unique(c(do_results_m$CellType,do_results_h$CellType))
Ncell = length(cell)

nodes <- data.frame( id=c(pheno_m,pheno_h,cell),
                     label=c(pheno_m,pheno_h,cell),
                     shape=c(rep("box",Npheno_m),rep("box",Npheno_h),rep("ellipse",Ncell)),
                     color=c(rep("lightblue",Npheno_m),rep("lightgreen",Npheno_h),rep("orange",Ncell)))

edges <- data.frame(from = c(adjacency$from),to = c(adjacency$to),arrows = c("to"))
visNetwork(nodes, edges)
}}
plot_pheno_cell = function(disease){
  if(disease != ""){
    if(disease %in% signif_results_h$Disease){
      do_results_m = signif_results_m[signif_results_m$MONDO_label == disease,]
      do_results_h = signif_results_h[signif_results_h$Disease == disease,]
    }else{
      disease = unique(signif_results_m[signif_results_m$DO.Disease == disease,"MONDO_label"])
      do_results_m = signif_results_m[signif_results_m$MONDO_label == disease,]
      do_results_h = signif_results_h[signif_results_h$Disease == disease,]
    }
      do_signif_counts = data.frame()
      for (c in unique(do_results_m$CellType)) {
        n_signif = length(do_results_m[do_results_m$CellType == c, ]$q)
        do_signif_counts = rbind(do_signif_counts,
                                 data.frame("PO"="MPO",
                                            "CellType"=c,
                                            "n_signif"=n_signif))}
      for (c in unique(do_results_h$CellType)) {
        n_signif = length(do_results_h[do_results_h$CellType == c, ]$q)
        do_signif_counts = rbind(do_signif_counts,
                                 data.frame("PO"="HPO",
                                            "CellType"=c,
                                            "n_signif"=n_signif))}
      
      do_signif_counts$CellType = factor(do_signif_counts$CellType,levels = cell_order$x)

      do_plt <- ggplot(do_signif_counts, aes(x=CellType,y=n_signif,fill=PO)) +
        geom_bar(stat="identity", position="stack",color="white",width=0.8,size=0.25)+
        cowplot::theme_cowplot()+
        ylab("N phenotypes") +
        scale_y_continuous(expand=c(0,0), limits= c(0,max(do_signif_counts$n_signif)+2))+
        theme(axis.text.x = element_text(angle=90,hjust=1,vjust = 0.3),axis.title.x = element_blank(),
              plot.title = element_text(hjust = 0.5),legend.title =element_blank()) +
        #coord_flip() +
        ggtitle(disease)
      return (do_plt)
    }
  }

table_do_m = function(disease){
  if(disease != ""){
    do_results_m = signif_results_m[signif_results_m$DO.Disease == disease |signif_results_m$MONDO_label == disease,]
      return(do_results_m)}}

table_do_h = function(disease){
  if(disease != ""){
    if(disease %in% signif_results_h$Disease){
      do_results_h = signif_results_h[signif_results_h$Disease == disease,]
    }else{
      disease = unique(signif_results_m[signif_results_m$DO.Disease == disease,"MONDO_label"])
      do_results_h = signif_results_h[signif_results_h$Disease == disease,]
    }
    return(do_results_h)}}