flat2Cdf<-function(file,chipType,tags=NULL,rows=1050,cols=1050,verbose=10,xynames=c("X","Y"),
    gcol=5,ucol=6,splitn=4,col.class=c("integer","character")[c(1,1,1,2,2,2)],...) {
  split.quick<-
  function(r,ucol,splitn=3,verbose=TRUE) {
    rn3<-substr(r[,ucol],1,splitn)
    split.matrix<-split.data.frame
    rr<-split(r,factor(rn3))
    if (verbose) cat(" split into",length(rr),"initial chunks ...")
    rr<-unlist(lapply(rr,FUN=function(u) split(u,u[,ucol])),recursive=FALSE)
    if (verbose) cat(" unwrapped into",length(rr),"chunks ...")
    names(rr)<-substr(names(rr),splitn+2,nchar(rr))
    rr
  }

  if (verbose) cat("Reading TXT file ...")
  file<-read.table(file,header=TRUE,colClasses=col.class,stringsAsFactors=FALSE,comment.char="",...)
  if (verbose) cat(" Done.\n")

  if (verbose) cat("Splitting TXT file indices into units ...")
  #gxys<-split.quick(file,ucol,splitn)
  gxysInd<-split(seq_len(nrow(file)),file[,ucol])
  if (verbose) cat(" Done.\n")

  l<-vector("list",length(gxysInd))
  if (verbose) cat("Creating structure for",length(gxysInd),"units (dot=250):\n")
  for(i in  1:length(gxysInd)) {
    thisTab <- file[ gxysInd[[i]], ]
    sp<-split(thisTab, factor(thisTab[,gcol]))
    #sp<-split(gxys[[i]],factor(gxys[[i]][,gcol]))
    e<-vector("list",length(sp))
    for(j in 1:length(sp)) {
      np<-nrow(sp[[j]])
      e[[j]]<-list(x=sp[[j]][,xynames[1]],y=sp[[j]][,xynames[2]],pbase=rep("A",np),tbase=rep("T",np),atom=0:(np-1),indexpos=0:(np-1),
                   groupdirection="sense",natoms=np,ncellsperatom=1)
    }
    names(e)<-names(sp)
    #l[[i]]<-list(unittype=1,unitdirection=1,groups=e,natoms=nrow(gxys[[i]]),ncells=nrow(gxys[[i]]),ncellsperatom=1,unitnumber=i)
    l[[i]]<-list(unittype=1,unitdirection=1,groups=e,natoms=nrow(thisTab),ncells=nrow(thisTab),ncellsperatom=1,unitnumber=i)
    if (verbose) { if(i %% 250==0) cat("."); if(i %% 5000==0) cat("(",i,")\n",sep="") }
  }
  rm(file,e,sp,thisTab); gc()
  cat("\n")
  #names(l)<-names(gxys)
  names(l)<-names(gxysInd)  
  rm(gxysInd); gc()
  if(!is.null(tags) && tags!="") filename<-paste(chipType,tags,sep=".")
  else filename<-chipType
  filename<-paste(filename,"cdf",sep=".")
  hdr<-list(probesets=length(l),qcprobesets=0,reference="",chiptype=chipType,filename=filename,
            nqcunits=0,nunits=length(l),rows=rows,cols=cols,refseq="",nrows=rows,ncols=cols)
  require(affxparser)
  writeCdf(hdr$filename, cdfheader=hdr, cdf=l, cdfqc=NULL, overwrite=TRUE, verbose=verbose)
  invisible(list(cdfList=l,cdfHeader=hdr))
}