GMPR<-function(OTUmatrix, min_ct=2, intersect_no=4){
  #datatype check
  if(!(class(OTUmatrix) %in% c("data.frame","matrix")))
    stop("Unknown datatype of object \"OTUmatrix\".")
  rownames(OTUmatrix)->SampleName
  if(ncol(OTUmatrix)<nrow(OTUmatrix))
    warning("Sample size is larger than OTU number. Check if samples are arranged in columns.")
  if(length(OTUmatrix[OTUmatrix<1 & OTUmatrix>0])>=0.1*ncol(OTUmatrix)*nrow(OTUmatrix))
    stop("More than 10% values are fractional, please check.")
  apply(OTUmatrix, MARGIN = 2, as.integer)->OTUmatrix

  gmpr(OTUmatrix,min_ct,intersect_no)->size.factor
  names(size.factor)<-SampleName
  size.factor[abs(size.factor-1)<1e-10]<-NA
  size.factor
}
