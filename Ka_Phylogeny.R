##R code for phylogeny comparison of fasta sequences
library(ape)
library(seqinr)
library(msa)
library(ggtree)

X = "file_sequences.fasta"
X_strings= readAAStringSet(X)
X_ClustalW = msa(X_strings, method = "ClustalW")
X_ClustalW2 = msaConvert(X_ClustalW, type = "seqinr::alignment")
XW <- dist.alignment(X_ClustalW2, "identity")
X_Tree = nj(XW)

msaplot(ggtree(X_Tree, branch.length = 'none',ndigits=3) + geom_text2(aes(subset=!isTip, label = ""), hjust=-.03) +
          geom_text(aes(x=branch, label= round(node.depth.edgelength(X_Tree),3)), vjust=0, color="firebrick") + 
          geom_tiplab(), fasta = X ,  offset=8)

msaPrettyPrint(X_ClustalW, output=c("pdf", "tex", "dvi", "asis"),  subset=NULL, file=NULL, alFile=NULL,
               askForOverwrite=TRUE,  psFonts=FALSE, code=NA,
               paperWidth=11, paperHeight=8.5, margins=c(0.1, 0.3),
               shadingMode="identical", shadingModeArg=NA, shadingColors="blues",showConsensus="bottom",
               consensusColors="ColdHot", consensusThreshold=50,showLogo="top", 
               logoColors="chemical",showLogoScale="none", showNames="left",
               showNumbering="right", showLegend=TRUE, furtherCode=NA, verbose=FALSE)
