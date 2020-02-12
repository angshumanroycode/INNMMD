innmmd.test.preset<-function(n,d,parameters=NULL,B.parameters=1000){
  if(is.null(parameters))
    parameters=PARAMETERS(n,d,B.parameters)
  return(INNMMDPRESET(n,d,parameters))
}

innmmd.test<-function(Xlist=NULL,Dlist=NULL,alpha=0.05,dist.type="euclidean",
                      lp=2,parameters=NULL,B.parameters=1000,B=100,
                      test.preset=NULL){
  if(!is.null(Xlist))
    Dlist=lapply(Xlist,function(x) dist(x,method=dist.type,p=lp))
  Dlist=lapply(Dlist,as.matrix)
  ns=as.numeric(sapply(Dlist,dim))
  n=ns[1]
  if(sum(n!=ns)) stop("Invalid Input")
  d=length(Dlist)
  if(is.null(test.preset))
    test.preset=innmmd.test.preset(n,d,parameters,B.parameters)
  result0=INNMMDTEST(Dlist,n,d,test.preset,B)
  val=result0[[1]]
  permval=result0[[2]]
  cutoff=apply(permval,1,function(x) as.numeric(quantile(x,1-alpha)))
  pval=(rowSums(permval>val)+1)/(B+1)
  result=list(val[1],val[2],cutoff[1],cutoff[2],pval[1],pval[2])
  names(result)=c("Tsum.stat","Tmax.stat","Tsum.cutoff","Tmax.cutoff",
                  "Tsum.pvalue","Tmax.pvalue")
  return(result)
}