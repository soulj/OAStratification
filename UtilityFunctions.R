#Utility functions

pSVA <- function(D, batches,num) {
  batch.mod <- model.matrix(~factor(batches))
  rSVA.SV <- sva(dat=D, mod=batch.mod, n.sv=num)
  rSVA.fit <- lmFit(D, cbind(batch.mod,rSVA.SV$sv))
  rSVA.D <- sweep(rSVA.fit$coefficients[,colnames(rSVA.fit$coefficients)==""]%*%
                    t(rSVA.SV$sv),1,rSVA.fit$coefficients[,"(Intercept)"],FUN="+")
  colnames(rSVA.D) <- colnames(D)
  return(rSVA.D)
}

GTFtoGeneLength <- function (gtffile){
  suppressPackageStartupMessages(library("EnsDb.Hsapiens.v79"))
  suppressPackageStartupMessages(library("EnsDb.Mmusculus.v79"))
  suppressPackageStartupMessages(library("GenomicFeatures"))
  
  
  txdb <- makeTxDbFromGFF(gtffile, format = "gtf", 
                          circ_seqs = character() )
  ebg <- exonsBy(txdb, by = "gene")
  exonic.gene.sizes <- lapply(ebg, function(x){sum(width(reduce(x)))})
  exonic.gene.sizes <-stack(exonic.gene.sizes)
  colnames(exonic.gene.sizes)<-c("length","gene_id")
  return(exonic.gene.sizes)
  
}

getEnrichedReactomePathways<-function(data,sigData,geneLengths,threshold=0.05){
  
  data<-as.data.frame(data)
  sigData<-as.data.frame(sigData)
  
  genes<-ifelse(as.data.frame(data)[,1]%in% as.data.frame(sigData)[,1],1,0)
  names(genes)<-data[,1]
  
  #map from ens to REACTOME
  ens2eg<-as.list(org.Hs.egENSEMBL2EG)

  eg2reactome=as.list(reactomeEXTID2PATHID)
  grepREACTOME=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
  reactome=lapply(ens2eg,grepREACTOME,eg2reactome)
  
  geneLengths<-geneLengths[match(data[,1],geneLengths$gene_id),]
  
  x<-nullp(genes,bias.data = geneLengths$length,plot.fit=T)
  
  REACTOME=goseq(x,gene2cat=reactome)
  REACTOME$padj=p.adjust(REACTOME$over_represented_pvalue,method="BH")
  xx <- as.list(reactomePATHID2NAME)
  REACTOME$Term=apply(REACTOME,1,function(x) xx[[unlist(x[1])]])
  REACTOME$Enrichment=REACTOME$numDEInCat/REACTOME$numInCat*100
  REACTOME$Adjpvaluelog=-log10(REACTOME$padj)
  REACTOME$Term<-unlist(REACTOME$Term)
  REACTOME.sig=REACTOME[REACTOME$padj<=threshold,]
  
  if (nrow(REACTOME.sig)==0) {
    print("no significant pathways")
    return()
  }
  
  reactomeResults=list()
  
  for ( i in 1:nrow(REACTOME.sig)) {
    
    #search reactome for the reactome term of interest and filter by differentially expressed genes
    reactomeTerm=REACTOME.sig$category[i]
    index=sapply(reactome,function(x) reactomeTerm %in% x)
    termIDs=names(index[index=="TRUE"])
    
    sig=sigData[sigData[,1] %in% termIDs ,]
    reactomeResults[[reactomeTerm]]=sig[,"gene_name"]
  }
  names(reactomeResults)=REACTOME.sig$Term
  
  reactomeResults=lapply(reactomeResults,function(x) paste(x, sep="", collapse=" ") )
  reactomeResults=data.frame(Term= names(reactomeResults),Genes = unlist(reactomeResults),Adj.pvalue=REACTOME.sig$padj,Enrichment=REACTOME.sig$Enrichment,Adjpvaluelog=REACTOME.sig$Adjpvaluelog)
  reactomeResults[,1]<-gsub( "Homo sapiens: ", "", as.character(reactomeResults[,1]))
  
  return(reactomeResults)
  
}