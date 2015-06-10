featMat = function(bam,chrstats,mappability_file="wgEncodeCrgMapabilityAlign100mer.bw"){
 nm = strsplit(names(bam)[1],split=":")[[1]]
 chr = nm[1]
 st = as.numeric(strsplit(nm[2],split="-")[[1]][1])
 ed = as.numeric(strsplit(nm[2],split="-")[[1]][2])
 
 ### set up the windows we'll use throughout
 sq = seq(st,ed,100)
 ir = IRanges(start=sq,width=100)

 ### get the flags of interest
 fl = expFlags(bam[[1]]$flag)
 samestr = fl$strand == fl$mstrand
 qcfail = fl$qcfail | fl$dupl

 ### pull out pos, mpos, isize, and mrnm
 pos = bam[[1]]$pos[!qcfail]
 mpos = bam[[1]]$mpos[!qcfail]
 isize = bam[[1]]$isize[!qcfail]
 mrnm = bam[[1]]$mrnm[!qcfail]
 samestr = samestr[!qcfail]
 qw = bam[[1]]$qwidth[!qcfail]
 qn = bam[[1]]$qname[!qcfail]
 s = bam[[1]]$seq[!qcfail]
 p1 = pmin(pos,mpos)
 p2 = pmax(pos,mpos)

 ### index the pairs that map exclusively to 
 ### this chromosome
 thischr = mrnm==chr & !fl$munmapped[!qcfail] 
 thischr[is.na(thischr)] = FALSE

 ###################################################
 ### get coverage weighted by discordant read pair status;
 ### read pairs are binned by their (apparent) insert size
 ###################################################
 #ol = log10(median(abs(isize)[thischr],na.rm=T) + 3*mad(abs(isize)[thischr],na.rm=T))
 #browser()
 id = bam[[1]]$tag$RG
 if(length(id) < 1 & grepl(":",qn[1]) & length(chrstats$median)>1){
  id = strsplit(qn,split=":")
  id1 = sapply(id,'[',1)
  id2 = sapply(id,'[',2)
  id = as.factor(as.character(as.factor(id1):as.factor(id2)))
 }
 id = id[thischr]

 ## a fix for simulated data with no RG info (see overrideRG in chrStats)
 if(is.null(id)){ 
  id = as.factor(rep(1,sum(thischr)))
  ol = log10(chrstats$ol[rep(1,sum(thischr))])
 }else{
  ol = log10(chrstats$ol[id])
 }
 ##

 bin = cut(log10(abs(isize)[thischr]),breaks=seq(2,6.5,length.out=8))
 ol_weight = as.numeric(abs(isize)[thischr] > 10^ol & !samestr[thischr])

 ## wrap up the binned read pairs into IRanges objects
 ol_cvg = tapply(1:sum(thischr),bin,
   function(i) list(cvg=IRanges(start=p1[thischr][i],end=p2[thischr][i]+qw[thischr][i]),
	w=ol_weight[i]))

 ## this function gets the coverage
 f_cvg = function(x,width){
  if(is.null(x)) return(Rle(0,width))
  if(all(is.na(x$w))) return(Rle(0,width)) 
  coverage(x$cvg,width=width,weight=x$w)
 }
 ## ...apply it over the list
 ol_cvg = lapply(ol_cvg,f_cvg,width=ed)

 ## chop it up and get windowed coverage
 v = lapply(ol_cvg,function(x) Views(x,ir))
 vm = lapply(v,viewMeans,na.rm=T)
 ol_cvg = do.call('cbind',vm)
 
 ol_cvg = log2((ol_cvg+runif(length(ol_cvg),0.01,0.1))/chrstats[["median"]])


 #######################################
 ### get inversion coverage while we're at it
 #######################################
 #browser()
 ol_weight = as.numeric(abs(isize)[thischr] > 10^ol & samestr[thischr] )
 # inv_cvg = coverage(IRanges(start=p1[thischr],end=p2[thischr]+qw[thischr]),
 #	weight=ol_weight*samestr[thischr])
 ## wrap up the binned read pairs into IRanges objects
 inv_cvg = tapply(1:sum(thischr),bin,
   function(i) list(cvg=IRanges(start=p1[thischr][i],end=p2[thischr][i]+qw[thischr][i]),
	w=ol_weight[i]))
 inv_cvg = lapply(inv_cvg,f_cvg,width=ed)
 v = lapply(inv_cvg,function(x) Views(x,ir))
 vm = lapply(v,viewMeans,na.rm=T)
 inv_cvg = do.call('cbind',vm)
 
 inv_cvg = log2((inv_cvg+runif(length(inv_cvg),0.01,0.1))/chrstats[["median"]])



# inv_cvg = viewMeans(Views(inv_cvg,ir),na.rm=T)
# inv_cvg = log2((inv_cvg+runif(length(inv_cvg),0.01,0.1))/chrstats[["median"]])
 #######################################
 ### vanilla coverage at various scales
 #######################################
 cvg = coverage(IRanges(start=pos[thischr],width=qw[thischr]))
 v1 = Views(cvg,IRanges(start=seq(st,ed,100),width=100))
 v2 = Views(cvg,IRanges(start=seq(st,ed,1000),width=1000))
 v3 = Views(cvg,IRanges(start=seq(st,ed,10000),width=10000))
 v4 = Views(cvg,IRanges(start=seq(st,ed,100000),width=100000))
 v5 = Views(cvg,IRanges(start=seq(st,ed,1000000),width=1000000))

 ir2 = IRanges(start(ir)-st,end(ir)-st)
 v1 = viewMeans(Views(Rle(viewMeans(v1,na.rm=T),width(v1)),ir2),na.rm=T)
 v2 = viewMeans(Views(Rle(viewMeans(v2,na.rm=T),width(v2)),ir2),na.rm=T)
 v3 = viewMeans(Views(Rle(viewMeans(v3,na.rm=T),width(v3)),ir2),na.rm=T)
 v4 = viewMeans(Views(Rle(viewMeans(v4,na.rm=T),width(v4)),ir2),na.rm=T)
 v5 = viewMeans(Views(Rle(viewMeans(v5,na.rm=T),width(v5)),ir2),na.rm=T)
 #print("chkpt 4")
 #browser()
 cvg = cbind(v1,v2,v3,v4,v5)
 cvg[is.na(cvg)] = 0

 #########################################
 ### GC-content coverage
 #########################################
 w = alphabetFrequency(s[thischr])
 w = (w[,2]+w[,3])/qw[thischr]

 gccvg = coverage(IRanges(start=pos[thischr],width=qw[thischr]),weight=w)
 gccvg = viewMeans(Views(gccvg,ir),na.rm=T)
 gccvg = gccvg/(v1+1)
 gccvg = log2(gccvg + runif(length(gccvg),0.01,0.1))
 


 ### normalize *vanilla* coverage & put on log2 scale
 cvg = log2((cvg+runif(length(cvg),0.01,0.1))/chrstats[["median"]])


 #########################################
 ### add the mappability scores
 #########################################
 require(rtracklayer)
 #browser()
 if(grepl("chr",chr)[1]){
  rg = RangedData(IRanges(start=st,end=ed),space=chr)
  rg = as(rg,"GRanges")
  mp = import(mappability_file, which = rg)
  mp_cvg = coverage(mp,weight="score",width=ed)[[chr]]
 }else{
  rg = RangedData(IRanges(start=st,end=ed),space=paste("chr",chr,sep=""))
  rg = as(rg,"GRanges")
  mp = import(mappability_file, which = rg)
  mp_cvg = coverage(mp,weight="score",width=ed)[[paste("chr",chr,sep="")]]
 }
 mp_cvg = viewMeans(Views(mp_cvg,ir))
 #print("chkpt 5")
 #browser()
 X = cbind(ol_cvg,inv_cvg,cvg,gccvg,mp_cvg)
 colnames(X) = c(paste("ol_cvg_",c(1:7),sep=""),paste("inv_cvg_",c(1:7),sep=""),
	paste("cvg_",c(1:5),sep=""),"GC_content",
	"mappability")
 
 out = RangedData(ir,space=chr,X)
 return(out)

}

###############################################################
### chrStats
###############################################################
chrStats = function (file, chr,recreateRG=FALSE,verbose=TRUE,overrideRG=FALSE) {
    if (length(chr) > 1) 
        stop("only meant for 1 chromosome at a time")
    h = scanBamHeader(file)[[1]]
    chrlength = h$targets[chr]
    if(chrlength > 5e4){
     st = seq(1, chrlength, length.out=200)
     ed = st + 10000
    }else{
     st = 1
     ed = chrlength
    }
    ir = RangesList(IRanges(start = st, end = ed))
    names(ir)[1] = chr
    
    dat = getBAM(file, ir, what = c("rname", "qname", "pos", 
        "qwidth", "isize"))
    rg = unlist(sapply(dat,function(x) x$tag$RG))

    ## check for RG; if empty, stop unless recreateRG is TRUE
    if(length(rg)>0 | recreateRG){
     if(recreateRG & length(rg)<1){ 
	if(verbose) warning("attempting to reconstruct lane assignments (pseudo-RG)")
	id = strsplit(unlist(sapply(dat,'[[',"qname")),split=":")
	id1 = sapply(id,'[',1)
	id2 = sapply(id,'[',2)
	id = as.factor(as.character(as.factor(id1):as.factor(id2)))
     }else{
      id = rg
     }     
    }else{
      if(!overrideRG) stop("no read group information in the BAM file; you may attempt to recreate it (EXPERIMENTAL) by setting recreateRG to TRUE")

    }
    #browser()
    pos = unlist(sapply(dat, "[[", "pos"))
    wd = unlist(sapply(dat, "[[", "qwidth"))
    covg = coverage(IRanges(start = pos, width = wd))
    covg = Views(covg, ir[[1]])
    m = do.call("c", viewApply(covg, "["))
    m = median(m[m>0])

    isize = abs(unlist(sapply(dat, "[[", "isize")))
    if(overrideRG) id = as.factor(rep(1,length(isize)))
    ol = tapply(abs(isize),id,median,na.rm=T) 
    m_ol = 4*tapply(abs(isize),id,mad,na.rm=T)
    ol = ol+m_ol
    if(any(m_ol/4 > 100) & verbose) warning("your fragment size distribution seems quite wide; are you sure this is a single library? Your ability to detect smaller SVs may be compromised.")
    if(nlevels(id)>15 & verbose) warning("there seem to be quite a few libraries and/or lanes -- are you sure this a single sample?  Was it multiplexed with other samples on a single lane?")
    out = list(median=round(m,0),ol=round(ol,0)) 
    return(out)
}

##############################################################
### getBAM
##############################################################
getBAM = function (file, rg, what = c("qname", "pos", "qwidth", "isize", 
    "mpos", "strand", "mapq", "mrnm", "flag", "seq", "cigar"), 
    flag = scanBamFlag(isUnmappedQuery = FALSE)) 
{
    require(Rsamtools)
    param = ScanBamParam(which = rg, what = what, flag = flag, tag="RG")
    out = scanBam(file = file, index = file, param = param)
    return(out)
}


################################################################
### a wrapper to do the segmentation
################################################################
seg = function(x,rf,rf_CN,fit_CN){
 require(HiddenMarkov)
 require(forestSV)
 load(x)
 tm = matrix(c(1-1e-9,1e-9,1e-9,1-1e-9),2,2)
 #pr = c(2.77, 459.49, 6.48, 4.68) ## old parameters
 pr_del = c(4.3111243,134.4780046,5.0147405,0.5178784)
 pr_dup = c(1.7533146,35.3081071,10.5557207,0.7526551)
 pr_inv = c(29.8835061,1260.5957387,2.4338932,0.4161289)
 nd = as.matrix(na.roughfix(as.data.frame(X)[,-c(1:4)]))
 prd = predict(rf,nd,type='prob')
 prd = prd + runif(length(prd),0.00001,0.01)
 prd[prd>0.99999] = 0.99999
 cn = exp(predict(fit_CN,predict(rf_CN,cbind(nd,prd),type='prob')))
 hmm_del = dthmm(prd[,1],
	tm,c(0.5,0.5),"beta",pm=list(shape1=c(pr_del[1],pr_del[3]),shape2=c(pr_del[2],pr_del[4])))
 st_del = Viterbi(hmm_del)
 hmm_dup = dthmm(prd[,2],
	tm,c(0.5,0.5),"beta",pm=list(shape1=c(pr_dup[1],pr_dup[3]),shape2=c(pr_dup[2],pr_dup[4])))
 st_dup = Viterbi(hmm_dup)
 hmm_inv = dthmm(prd[,4],
	tm,c(0.5,0.5),"beta",pm=list(shape1=c(pr_inv[1],pr_inv[3]),shape2=c(pr_inv[2],pr_inv[4])))
 st_inv = Viterbi(hmm_inv)

 out = RangedData(ranges(X)[[1]],space=space(X)[1],
	data.frame(prd,CN=cn,del_state=st_del==2,dup_state=st_dup==2,inv_state=st_inv==2))
 return(out)
}

### a wrapper to merge the calls and get mean scores/CN per call
seg2 = function(x){
 x = do.call('rbind',x)
 chr = unique(space(x))
 if(length(chr)>1) stop("can only work on one chromosome at a time")
 oo = order(start(x))
 x = x[oo,] ## sort the RD object
 idx = seq(1,max(end(x)),100) ## make sure we have all starts

 ### get the coords of the calls
 del_rg = ranges(reduce(x[x$del_state,]))[[1]]
 dup_rg = ranges(reduce(x[x$dup_state,]))[[1]]
 inv_rg = ranges(reduce(x[x$inv_state,]))[[1]]
 cpx_rg = ranges(reduce(x[x$del_state & x$dup_state,]))[[1]]

 ### create chromosome-length score vectors
 del_s = numeric(length(idx))
 del_s[match(start(x),idx)] = x$deletion
 del_s = Rle(del_s,rep(100,length(idx)))

 dup_s = numeric(length(idx))
 dup_s[match(start(x),idx)] = x$duplication
 dup_s = Rle(dup_s,rep(100,length(idx)))

 inv_s = numeric(length(idx))
 inv_s[match(start(x),idx)] = x$inversion
 inv_s = Rle(inv_s,rep(100,length(idx)))

 cn_s = numeric(length(idx))
 cn_s[match(start(x),idx)] = x$CN
 cn_s = Rle(cn_s,rep(100,length(idx)))

 ### get the average score and CN per call
 delsS = viewMeans(Views(del_s,del_rg))
 delsCN = viewMeans(Views(cn_s,del_rg))
 cpx = del_rg%over%cpx_rg
 DEL = RangedData(del_rg,space=rep(chr,length(del_rg)),
		score=delsS,CN=delsCN,type=rep("DEL",length(del_rg)),complex=cpx)

 dupsS = viewMeans(Views(dup_s,dup_rg))
 dupsCN = viewMeans(Views(cn_s,dup_rg))
 cpx = dup_rg%over%cpx_rg
 DUP = RangedData(dup_rg,space=rep(chr,length(dup_rg)),
	score=dupsS,CN=dupsCN,type=rep("DUP",length(dup_rg)),complex=cpx)

 invS = viewMeans(Views(inv_s,inv_rg))
 invCN = viewMeans(Views(cn_s,inv_rg))
 cpx = inv_rg%over%cpx_rg
 INV = RangedData(inv_rg,space=rep(chr,length(inv_rg)),
	score=invS,CN=invCN,type=rep("INV",length(inv_rg)),complex=cpx)

 return(rbind(DEL,DUP,INV))

}

### the outermost wrapper
segSV = function(dir,rf,rf_CN,fit_CN,cl){
 ff = list.files(dir,pattern="features_",full.names=T)
 chr = unique(sapply(strsplit(ff,"_"),'[',2))
 if(length(chr)<1) stop("where are the feature files?")
 ll = list()
 for(i in 1:length(chr)){
  ff = list.files(dir,pattern=paste("features_",chr[i],"_",sep=""),full.names=T)
  if(length(ff)>1){
   tmp = parLapply(cl,ff,seg,rf,rf_CN,fit_CN)
  }else{
   tmp = list(seg(ff,rf,rf_CN,fit_CN))
  }
  ll[[i]] = seg2(tmp)
  print(paste("done with chromosome",chr[i]))
 }
 ll = do.call('rbind',ll)
 return(ll)
}

############################################################
### function to compute and save all features, given the path to a BAM file
############################################################
allFeatures = function(fn,basename,
	mf="/media/jmichaelson/scratch2/forestSV/variants/wgEncodeCrgMapabilityAlign100mer.bw",
	recreateRG=FALSE,overrideRG=FALSE){
 sl = scanBamHeader(fn)[[1]]$targets
 ## for now we have to limit to these chromosomes because these are the only ones 
 ## in the mappability file
 sl = sl[names(sl)%in%c(paste("chr",c(1:22,"X","Y","M"),sep=""),c(1:22,"X","Y"))]

 ll = lapply(sl,breakInChunks,1e7)
 chr = rep(names(ll),times=sapply(ll,length))
 ll = do.call('rbind',lapply(ll,as.data.frame))
 ll = data.frame(space=chr,ll,row.names=NULL)
 ll = as(ll,"RangedData")

 chr_ll = list()
 chr = names(sl)
 for(i in 1:length(chr)){
  chr_ll[[i]] = chrStats(fn,chr[i],recreateRG=recreateRG,overrideRG=overrideRG)
  print(paste(chr_ll[[i]]$median,chr[i]))
  gc()
 }
 names(chr_ll) = chr

 system(paste("mkdir",basename))
 for(i in 1:nrow(ll)){
  #mf = "/media/scratch2/forestSV/variants/wgEncodeCrgMapabilityAlign100mer.bw"
  rg = ll[i,]
  chr = as.character(space(rg)[1])
  bam = getBAM(fn,rg)
  if(length(bam[[1]]$pos)<10) next
  chrstats = chr_ll[[chr]]
  X = featMat(bam,chrstats,mf)
  save(X,file=paste(basename,"/features_",chr,"_",start(rg)[1],"_",end(rg)[1],".Rdata",sep=""))
  print(i)
  gc()
 }

}


#####################################################################
### reciprocal overlap
#####################################################################
rol = function(a,b){
 require(IRanges)
 ol = IRanges::findOverlaps(a,b)
 DF = data.frame()
 for(i in 1:length(ol)){
  ind = cbind(IRanges::queryHits(ol[[i]]),IRanges::subjectHits(ol[[i]]))
  A = a[IRanges::space(a)==names(ol)[i],]
  B = b[IRanges::space(b)==names(ol)[i],]
  if(nrow(ind) < 1) next
  
  so = pmin(IRanges::start(A[ind[,1],]),IRanges::start(B[ind[,2],]))
  si = pmax(IRanges::start(A[ind[,1],]),IRanges::start(B[ind[,2],]))
  eo = pmax(IRanges::end(A[ind[,1],]),IRanges::end(B[ind[,2],]))
  ei = pmin(IRanges::end(A[ind[,1],]),IRanges::end(B[ind[,2],]))
  S = (ei-si)/(eo-so)
  df = data.frame(chr=names(ol)[i],startA=IRanges::start(A[ind[,1],]),endA=IRanges::end(A[ind[,1],]),
        startB=IRanges::start(B[ind[,2],]),endB=IRanges::end(B[ind[,2],]),overlap=S)
  if("type"%in%c(colnames(A),colnames(B))){

   df$type=paste(A$type[ind[,1]],B$type[ind[,2]],sep="|")
  }
  df$scoreA=IRanges::score(A)[ind[,1]]
  df$scoreB=IRanges::score(B)[ind[,2]]
  DF = rbind(DF,df)

 }
 return(DF)
}

