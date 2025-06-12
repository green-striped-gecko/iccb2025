#### 

library(ShortRead )

ref =   readFasta("d:/bernd/r/Elise_pansnp/final.genome.scf.fasta")

#gl needs to have run via blast as it need sstart, send and sacc

gl.extract.sequence <- function(x, genome=ref, leftextend=0, rightextend=0, file.name=NULL, file.path=tempdir()) {
  
  saveseq <- DNAStringSet()
  for ( i in 1:nLoc(x)) {
  seqstart <- x@other$loc.metrics$sstart[i]
  seqend <- x@other$loc.metrics$send[i]
  rev <- FALSE
  if (seqstart > seqend) {
   dummy <- seqstart
   seqstart <- seqend
   seqend <- dummy
   dummy <- leftextend
   leftextend <- rightextend
   rightextend <- dummy
   rev <- TRUE
  }
  
  index <- which(as.character(id(ref))==x@other$loc.metrics$sacc[i])
  exseq <- narrow(ref[index], start = seqstart-leftextend, end = seqend+rightextend)
  if (rev)  exseq <- reverseComplement(sread(exseq)) else exseq <- sread(exseq)
  saveseq[i] <- exseq
  x@other$loc.metrics$exseq[i] <- as.character(exseq)
  x@other$loc.metrics$exseqstart[i] <- seqstart-leftextend
  x@other$loc.metrics$exseqend[i] <- seqend+rightextend
  }
  if (!is.null(file.name)) 
  {
    names(saveseq) <- x@other$loc.metrics$sacc
    writeFasta(saveseq, file=file.path(file.path,file.name))
  }
  return(x)
  
}

#function to extract sequences.... (at the moment the reverseComplement is shown)
y <-gl.extract.sequence(x, genome=ref, leftextend=10, rightextend=10, file.name="xyz.fasta", file.path="d:/temp")
(DNAStringSet(y@other$loc.metrics$TrimmedSequence))
(DNAStringSet(y@other$loc.metrics$exseq))
  
  