options(stringsAsFactors = FALSE, digits = 16, width = 180)
suppressPackageStartupMessages({library(edgeR); library(limma); library(tximport); library(jsonlite)})
for (p in c("limma", "edgeR", "tximport", "statmod", "jsonlite")) {
  stopifnot(requireNamespace(p, quietly = TRUE))
}
Sys.setlocale("LC_COLLATE", "C")
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
manifest <- read.delim("data/input_manifest.tsv", check.names = FALSE)
manifest <- manifest[startsWith(manifest$file, "data/airway/"), ]
stopifnot(nrow(manifest) > 0, identical(unname(tools::md5sum(manifest$file)), manifest$md5))

analyse_culture <- function() {
 options(stringsAsFactors=FALSE,digits=15,width=160)
 Sys.setlocale('LC_COLLATE','C'); Sys.setenv(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
 ROOT <- 'data/airway/culture'
 OUT <- 'output/tables/airway/culture'
 dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
 rt<-function(p)read.delim(p,check.names=FALSE,comment.char='',stringsAsFactors=FALSE)
 wt<-function(x,n)write.table(x,file.path(OUT,n),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
 wm<-function(x,n,id='gene_id'){
  d<-data.frame(id=rownames(x),x,check.names=FALSE);names(d)[1]<-id
  con<-gzfile(file.path(OUT,n),'wt');on.exit(close(con));write.table(d,con,sep='\t',row.names=FALSE,quote=TRUE,na='NA')
 }
 log<-function(s){s<-paste(format(Sys.time(),tz='UTC',usetz=TRUE),s);cat(s,'\n');write(s,file.path(OUT,'run_log.txt'),append=TRUE)}
 warns<-character();context<-'setup';status<-'FAILED'
 tryCatch(withCallingHandlers({
  for(p in c('limma','edgeR','tximport','statmod','jsonlite'))stopifnot(requireNamespace(p,quietly=TRUE))
  RNGkind('Mersenne-Twister','Inversion','Rejection')
  IN <- ROOT
  sm<-rt(file.path(IN,'sample_manifest.tsv'))
  libs<-c(sprintf('S%02d',1:6),sprintf('S%02d',13:18));donors<-sprintf('P%02d',1:6)
  stopifnot(identical(sm$library,libs),all(table(sm$donor,sm$condition)==1),
            identical(sm$donor,rep(donors,each=2)),identical(sm$group,rep(c('PC','LC'),each=6)),
            identical(sm$condition,rep(c('baseline','stimulation'),6)))
  qc <- rt(file.path(ROOT, 'technical_QC.tsv')); stopifnot(all(qc$status == 'TECHNICAL_PASS'))
  wt(qc, 'technical_QC.tsv'); wt(sm, 'sample_metadata.tsv')
  design<-function(s){
   d<-factor(s$donor);X<-model.matrix(~0+d)
   colnames(X)<-levels(d)
   X<-cbind(X,PC_stim=as.numeric(s$group=='PC'&s$condition=='stimulation'),
     LC_stim=as.numeric(s$group=='LC'&s$condition=='stimulation'))
   rownames(X)<-s$library;stopifnot(qr(X)$rank==ncol(X));X
  }
  X<-design(sm);ct<-setNames(c(rep(0,6),-1,1),colnames(X))
  stopifnot(ncol(X)==8,nrow(X)-qr(X)$rank==4)
  wt(data.frame(library=rownames(X),X,check.names=FALSE),'primary_design.tsv')
  wt(data.frame(coefficient=names(ct),contrast=ct),'primary_contrast.tsv')
  context<-'import';log('Importing 12 paired culture libraries.')
  tx2g<-rt(file.path(IN,'tx2gene.tsv'))
  stopifnot(nrow(tx2g)==252835,!anyDuplicated(tx2g$transcript_id),length(unique(tx2g$gene_id))==62700)
  files<-setNames(file.path(IN,'quants',libs,'quant.sf.gz'),libs)
  for(f in files){q<-rt(gzfile(f));stopifnot(identical(q$Name,tx2g$transcript_id),
    all(is.finite(as.matrix(q[,-1]))),all(as.matrix(q[,-1])>=0),abs(sum(q$TPM)-1e6)<2)}
  txi<-tximport::tximport(files,type='salmon',tx2gene=tx2g[,1:2],countsFromAbundance='lengthScaledTPM')
  stopifnot(identical(colnames(txi$counts),libs),nrow(txi$counts)==62700,
    all(is.finite(txi$counts)),all(txi$counts>=0),txi$countsFromAbundance=='lengthScaledTPM')
  saveRDS(txi,file.path(OUT,'tximport_lengthScaledTPM.rds'),compress='xz')
  g<-tx2g[!duplicated(tx2g$gene_id),c('gene_id','gene_symbol_gencode')]
  g<-g[match(rownames(txi$counts),g$gene_id),];stopifnot(identical(g$gene_id,rownames(txi$counts)))
  hg<-rt(gzfile(file.path('data/airway/assets','hgnc_complete_set_20261017.txt.gz')))
  hg<-hg[hg$status=='Approved',];prev<-strsplit(ifelse(is.na(hg$prev_symbol),'',hg$prev_symbol),'|',fixed=TRUE)
  own<-split(rep(hg$symbol,lengths(prev)),unlist(prev,use.names=FALSE))
  prevmap<-vapply(own,function(v){u<-unique(v);if(length(u)==1)u else ''},character(1))
  canon<-function(s){
   out<-rep('',length(s));exact<-s%in%hg$symbol;out[exact]<-s[exact];j<-match(s[!exact],names(prevmap))
   out[!exact]<-ifelse(is.na(j),'',unname(prevmap[j]));out
  }
  g$symbol<-canon(g$gene_symbol_gencode)
  g$mapping_method<-ifelse(g$symbol==g$gene_symbol_gencode,'approved',ifelse(nzchar(g$symbol),'unique_previous','unresolved'))
  g$mean_TPM<-rowMeans(txi$abundance);g$mean_imported_count<-rowMeans(txi$counts)
  g$n_reference_ids_per_symbol<-ifelse(nzchar(g$symbol),as.integer(table(g$symbol)[g$symbol]),NA_integer_)
  y0<-edgeR::DGEList(txi$counts)
  keep<-edgeR::filterByExpr(y0,group=interaction(sm$group,sm$condition),min.count=10,min.total.count=15,large.n=10,min.prop=.7)
  g$primary_filter<-keep
  g$strict_filter<-keep&rowSums(edgeR::cpm(y0)>=1)>=6
  ord<-order(-g$mean_TPM,g$gene_id);ord<-ord[keep[ord]&nzchar(g$symbol[ord])]
  reps<-ord[!duplicated(g$symbol[ord])]
  g$set_representative<-seq_len(nrow(g))%in%reps
  wt(g,'gene_mapping_filter_details.tsv')
  panel<-rt(file.path('data/airway/assets','gene_set_panel.tsv'));members<-rt(file.path('data/airway/assets','gene_set_members.tsv'))
  members$symbol<-canon(members$gene_symbol)
  members$gene_id<-g$gene_id[reps][match(members$symbol,g$symbol[reps])]
  members$retained_primary<-!is.na(members$gene_id)&nzchar(members$symbol)
  wt(members,'gene_set_member_details.tsv')
  counts<-txi$counts[keep,,drop=FALSE]
  normalize_voom<-function(M,D,quality=FALSE){
   y<-edgeR::calcNormFactors(edgeR::DGEList(M),method='TMM')
   v<-if(quality)limma::voomWithQualityWeights(y,D,span=.5,normalize.method='none',plot=FALSE)
      else limma::voom(y,D,span=.5,normalize.method='none',plot=FALSE,save.plot=TRUE)
   stopifnot(identical(dim(v$weights),dim(v$E)))
   dimnames(v$weights)<-dimnames(v$E)
   list(y=y,v=v)
  }
  nv<-normalize_voom(counts,X);v<-nv$v
  wt(data.frame(library=libs,nv$y$samples,check.names=FALSE),'normalisation_factors.tsv')
  wm(v$E,'normalised_logCPM.tsv.gz');wm(v$weights,'voom_weights.tsv.gz')
  wm(txi$counts,'gene_counts_lengthScaledTPM.tsv.gz');wm(txi$abundance,'gene_TPM.tsv.gz')
  saveRDS(list(txi=txi,g=g,sm=sm,panel=panel,members=members,v=v,y=nv$y,X=X,ct=ct),file.path(OUT,'analysis_inputs.rds'),compress='xz')
  log(sprintf('Retained %d / %d gene IDs; %d unique mapped representatives.',sum(keep),length(keep),length(reps)))
  make_indices<-function(gids,uniqueonly=FALSE){
   setNames(lapply(panel$gene_set,function(gs){
    z<-unique(members$gene_id[members$gene_set==gs&members$retained_primary])
    if(uniqueonly)z<-intersect(z,g$gene_id[g$n_reference_ids_per_symbol==1&!is.na(g$n_reference_ids_per_symbol)])
    which(gids%in%z)
   }),panel$gene_set)
  }
  allrows<-list();alldonors<-list();modelmetadata<-list()
  run_test<-function(V,D,C,s,scenario,seed,uniqueonly=FALSE,genes=FALSE){
   context<<-scenario;log(paste('Fit',scenario));idx<-make_indices(rownames(V$E),uniqueonly)
   n<-lengths(idx);eligible<-n>=10&n/panel$set_size>=.2
   coverage<-data.frame(gene_set=panel$gene_set,source_size=panel$set_size,n_retained=n,coverage=n/panel$set_size,eligible=eligible)
   wt(coverage,paste0(scenario,'_coverage.tsv'))
   set.seed(seed)
   rr<-limma::mroast(V,index=idx[eligible],design=D,contrast=C,set.statistic='mean',nrot=19999L,
     robust=TRUE,trend=FALSE,approx.zscore=TRUE,legacy=FALSE,midp=FALSE,sort='none')
   rr<-data.frame(gene_set=rownames(rr),rr,check.names=FALSE)
   rr<-merge(coverage,rr,by='gene_set',all.x=TRUE,sort=FALSE);rr<-rr[match(panel$gene_set,rr$gene_set),]
   rr$FDR<-p.adjust(rr$PValue,method='BH',n=10);rr$FDR.Mixed<-p.adjust(rr$PValue.Mixed,method='BH',n=10)
   ca<-limma::contrastAsCoef(D,C)
   fit<-limma::eBayes(limma::lmFit(V,ca$design),robust=TRUE,trend=FALSE)
   ef<-fit$coefficients[,ca$coef]
   rr$mean_weighted_gene_logFC<-vapply(idx,function(j)if(length(j))mean(ef[j])else NA_real_,numeric(1))
   rr$median_weighted_gene_logFC<-vapply(idx,function(j)if(length(j))median(ef[j])else NA_real_,numeric(1))
   rr$scenario<-scenario;rr$seed<-seed;rr$residual_df<-nrow(D)-ncol(D)
   rr$n_donors<-length(unique(s$donor));rr$MonteCarlo_SE_approx<-sqrt(rr$PValue*(1-rr$PValue)/20000)
   rr$reason<-ifelse(eligible,'tested','insufficient gene coverage')
   wt(rr,paste0(scenario,'_ROAST.tsv'));allrows[[scenario]]<<-rr
   if(all(table(s$donor,s$condition)==1)&length(unique(s$condition))==2){
    pairs<-unique(s$donor);B<-match(paste(pairs,'baseline'),paste(s$donor,s$condition))
    S<-match(paste(pairs,'stimulation'),paste(s$donor,s$condition))
    de<-V$E[,S,drop=FALSE]-V$E[,B,drop=FALSE];colnames(de)<-pairs
    dd<-do.call(rbind,lapply(seq_along(idx),function(i)data.frame(scenario=scenario,gene_set=names(idx)[i],donor=pairs,group=s$group[B],
      mean_paired_logCPM=colMeans(de[idx[[i]],,drop=FALSE]),n_genes=length(idx[[i]]))))
    alldonors[[scenario]]<<-dd
    if(scenario=='primary'){wm(de,'gene_paired_logCPM_changes.tsv.gz');wt(dd,'primary_donor_program_changes.tsv')}
   }
   if(genes){
    gt<-limma::topTable(fit,coef=ca$coef,number=Inf,sort.by='none',confint=.95)
    gt<-data.frame(gene_id=rownames(gt),symbol=g$symbol[match(rownames(gt),g$gene_id)],gt,check.names=FALSE)
    wt(gt,paste0(scenario,'_gene_statistics.tsv'))
   }
   wt(data.frame(library=rownames(D),D,check.names=FALSE),paste0(scenario,'_design.tsv'))
   modelmetadata[[scenario]]<<-list(rank=qr(D)$rank,n_samples=nrow(D),residual_df=nrow(D)-ncol(D),contrast=C,
     eligible_sets=sum(eligible),seed=seed,nrot=19999,genes=nrow(V$E))
   if(scenario=='primary')saveRDS(list(v=V,design=D,contrast=C,fit=fit,contrast_coefficient=ca$coef,indices=idx),file.path(OUT,'primary_model.rds'),compress='xz')
   invisible(list(result=rr,indices=idx,fit=fit))
  }
  primary<-run_test(v,X,ct,sm,'primary',2026102601L,genes=TRUE)
  cp<-ct*0;cp['PC_stim']<-1;cl<-ct*0;cl['LC_stim']<-1
  run_test(v,X,cp,sm,'secondary_PC_stimulation',2026102602L,genes=TRUE)
  run_test(v,X,cl,sm,'secondary_LC_stimulation',2026102603L,genes=TRUE)
  b<-sm$condition=='baseline';sb<-sm[b,];db<-model.matrix(~0+factor(sb$group,levels=c('PC','LC')))
  colnames(db)<-c('PC','LC');rownames(db)<-sb$library;vb<-normalize_voom(counts[,b],db)$v
  run_test(vb,db,c(-1,1),sb,'secondary_baseline_LC_minus_PC',2026102604L,genes=TRUE)
  vq<-normalize_voom(counts,X,quality=TRUE)$v
  wt(data.frame(library=libs,vq$targets,check.names=FALSE),'quality_weight_sensitivity_samples.tsv')
  run_test(vq,X,ct,sm,'sensitivity_quality_weights',2026102611L)
  vh<-normalize_voom(txi$counts[g$strict_filter,],X)$v
  run_test(vh,X,ct,sm,'sensitivity_strict_filter',2026102612L)
  run_test(v,X,ct,sm,'sensitivity_unambiguous_symbols',2026102613L,uniqueonly=TRUE)
  for(i in seq_along(donors)){
   j<-sm$donor!=donors[i];ss<-sm[j,];D<-design(ss);C<-c(rep(0,5),-1,1)
   vv<-normalize_voom(counts[,j],D)$v
   run_test(vv,D,C,ss,paste0('omit_',donors[i]),2026102620L+i)
  }
  context<-'CAMERA';log('Competitive gene-set sensitivity analysis.')
  idx<-primary$indices;eligible<-lengths(idx)>=10&lengths(idx)/panel$set_size>=.2
  cam<-limma::camera(v,index=idx[eligible],design=X,contrast=ct,inter.gene.cor=NA,allow.neg.cor=FALSE,trend.var=FALSE,sort=FALSE)
  cam<-data.frame(gene_set=rownames(cam),cam);cam$FDR<-p.adjust(cam$PValue,'BH',n=10);wt(cam,'sensitivity_CAMERA.tsv')
  wt(do.call(rbind,allrows),'all_ROAST_results.tsv');wt(do.call(rbind,alldonors),'all_donor_program_changes.tsv')
  pc<-prcomp(t(v$E),center=TRUE,scale.=FALSE)
  wt(data.frame(sm,pc$x[,1:4]),'PCA_samples.tsv')
  wt(data.frame(PC=paste0('PC',seq_along(pc$sdev)),variance_fraction=pc$sdev^2/sum(pc$sdev^2)),'PCA_variance.tsv')
  wt(data.frame(gene_id=rownames(v$E),mean_logCPM=rowMeans(v$E),sd_logCPM=apply(v$E,1,sd),mean_weight=rowMeans(v$weights)), 'gene_expression_QC.tsv')
  saveRDS(list(modelmetadata=modelmetadata,voom_xy=v$voom.xy,voom_line=v$voom.line),file.path(OUT,'model_metadata.rds'))
  jsonlite::write_json(modelmetadata,file.path(OUT,'model_metadata.json'),auto_unbox=TRUE,pretty=TRUE)
  status<-'COMPLETE';log('All prespecified analyses complete.')
 },warning=function(w){warns<<-c(warns,paste(context,conditionMessage(w)));invokeRestart('muffleWarning')}),
 error=function(e){log(paste('FAILED',context,conditionMessage(e)));write(conditionMessage(e),file.path(OUT,'ERROR.txt'))})
 writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
 writeLines(unique(warns),file.path(OUT,'warnings.txt'))
 jsonlite::write_json(list(status=status,finished_utc=format(Sys.time(),tz='UTC',usetz=TRUE),context=context),file.path(OUT,'RUN_STATUS.json'),auto_unbox=TRUE,pretty=TRUE)
 if(status!='COMPLETE')stop('Analysis failed; review retained evidence.')
}

analyse_longitudinal <- function() {
options(stringsAsFactors=FALSE,digits=16,width=180)
suppressPackageStartupMessages({library(edgeR);library(limma);library(jsonlite)})
ROOT <- 'data/airway/longitudinal'; OUT <- 'output/tables/airway/longitudinal'
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
rt<-function(p)read.delim(p,check.names=FALSE,comment.char='')
wt<-function(x,n)write.table(x,file.path(OUT,n),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
wm<-function(x,n){con<-gzfile(file.path(OUT,n),'wt');write.table(data.frame(gene_id=rownames(x),x,check.names=FALSE),con,sep='\t',row.names=FALSE,quote=TRUE);close(con);z<-rt(gzfile(file.path(OUT,n)));stopifnot(nrow(z)==nrow(x),ncol(z)==ncol(x)+1)}
safe_save<-function(x,n){con<-gzfile(file.path(OUT,n),'wb',compression=6);saveRDS(x,con);close(con);stopifnot(identical(readRDS(file.path(OUT,n)),x))}
RNGkind('Mersenne-Twister','Inversion','Rejection')
warnings<-character();status<-'FAILED';error<-NULL
tryCatch(withCallingHandlers({
 meta<-rt(file.path(ROOT,'input/GSE261777_sample_metadata.tsv'))
 all<-read.delim(gzfile(file.path(ROOT,'input/GSE261777_counts.txt.gz')),dec=',',check.names=FALSE,row.names=1)
 all<-as.matrix(all);stopifnot(nrow(all)==34515,ncol(all)==55,all(is.finite(all)),all(all>=0),!anyDuplicated(rownames(all)))
 stopifnot(setequal(colnames(all),meta$sample_column),length(unique(meta$individual))==49)
 donors<-paste('patient',c(12,16,18,31,35,36));samples<-as.vector(rbind(c('12','16','18','31','35','36'),paste0(c('12','16','18','31','35','36'),'_2')))
 sm<-meta[match(samples,meta$sample_column),];sm$donor<-sm$individual;sm$late<-as.numeric(sm$visit=='study visit 2')
 stopifnot(all(table(sm$donor,sm$late)==1),identical(sm$donor,rep(donors,each=2)),all(sm$`disease state`=='long COVID'))
 wt(sm,'paired_sample_manifest.tsv');wt(meta,'all55_sample_manifest.tsv')
 C<-all[,samples,drop=FALSE]
 X<-model.matrix(~0+factor(sm$donor));colnames(X)<-donors;X<-cbind(X,late=sm$late);rownames(X)<-samples
 stopifnot(qr(X)$rank==7,nrow(X)-ncol(X)==5)
 contrast<-c(rep(0,6),1);names(contrast)<-colnames(X)
 wt(data.frame(sample=samples,X,check.names=FALSE),'primary_design.tsv')
 ref<-rt(file.path('data/airway/assets/gencode_v44_symbols.tsv'));ref$stable<-sub('[.][0-9]+$','',ref$gene_id)
 ref$symbol[is.na(ref$symbol)]<-''
 owners<-split(ref$symbol,ref$stable)
 canon<-vapply(owners,function(s){u<-unique(s[nzchar(s)]);if(length(u)==1)u else ''},character(1))
 g<-data.frame(gene_id=rownames(C),symbol=unname(canon[match(rownames(C),names(canon))]))
 g$symbol[is.na(g$symbol)]<-'';g$mapped<-nzchar(g$symbol)
 y0<-DGEList(C);keep<-filterByExpr(y0,group=sm$late,min.count=10,min.total.count=15,large.n=10,min.prop=.7)
 g$primary_filter<-keep;g$mean_CPM<-rowMeans(cpm(y0))
 ord<-order(-g$mean_CPM,g$gene_id);ord<-ord[keep[ord]&g$mapped[ord]];reps<-ord[!duplicated(g$symbol[ord])]
 g$set_representative<-seq_len(nrow(g))%in%reps;wt(g,'gene_mapping_filter_details.tsv')
 panel<-rt(file.path('data/airway/assets/gene_set_panel.tsv'));members<-rt(file.path('data/airway/assets/gene_set_members.tsv'))
 hg<-rt(gzfile(file.path('data/airway/assets/hgnc_complete_set_20261017.txt.gz')));hg<-hg[hg$status=='Approved',]
 pr<-strsplit(ifelse(is.na(hg$prev_symbol),'',hg$prev_symbol),'|',fixed=TRUE);own<-split(rep(hg$symbol,lengths(pr)),unlist(pr,use.names=FALSE))
 prev<-vapply(own,function(s){u<-unique(s);if(length(u)==1)u else ''},character(1))
 members$symbol<-ifelse(members$gene_symbol%in%hg$symbol,members$gene_symbol,unname(prev[match(members$gene_symbol,names(prev))]));members$symbol[is.na(members$symbol)]<-''
 members$gene_id<-g$gene_id[reps][match(members$symbol,g$symbol[reps])];members$retained<-!is.na(members$gene_id)
 wt(members,'gene_set_member_details.tsv')
 normalize<-function(c,d){y<-calcNormFactors(DGEList(c),method='TMM');v<-voom(y,d,span=.5,normalize.method='none',plot=FALSE,save.plot=TRUE);dimnames(v$weights)<-dimnames(v$E);list(y=y,v=v)}
 nv<-normalize(C[keep,,drop=FALSE],X);v<-nv$v;y<-nv$y
 wt(data.frame(sample=samples,y$samples,check.names=FALSE),'normalisation_factors.tsv')
 wm(C,'paired_submitted_counts.tsv.gz');wm(v$E,'normalised_logCPM.tsv.gz');wm(v$weights,'voom_weights.tsv.gz')
 idx<-setNames(lapply(panel$gene_set,function(gs)which(rownames(v$E)%in%members$gene_id[members$gene_set==gs&members$retained])),panel$gene_set)
 cov<-data.frame(gene_set=panel$gene_set,n_retained=lengths(idx),source_size=panel$set_size);cov$coverage<-cov$n_retained/cov$source_size;cov$eligible<-cov$n_retained>=10&cov$coverage>=.2;wt(cov,'program_coverage.tsv')
 run_roast<-function(V,D,CT,name,seed){
  set.seed(seed); rr<-mroast(V,index=idx[cov$eligible],design=D,contrast=CT,set.statistic='mean',nrot=19999L,robust=TRUE,trend=FALSE,approx.zscore=TRUE,legacy=FALSE,midp=FALSE,sort='none')
  rr<-data.frame(gene_set=rownames(rr),rr,check.names=FALSE);rr<-merge(cov,rr,by='gene_set',all.x=TRUE,sort=FALSE);rr<-rr[match(panel$gene_set,rr$gene_set),]
  rr$FDR<-p.adjust(rr$PValue,'BH',n=10);rr$FDR.Mixed<-p.adjust(rr$PValue.Mixed,'BH',n=10);rr$scenario<-name;rr$seed<-seed
  fit<-eBayes(lmFit(V,D),robust=TRUE,trend=FALSE);coef<-ncol(D)
  rr$mean_weighted_gene_logFC<-vapply(idx,function(i)mean(fit$coefficients[i,coef]),numeric(1));rr$n_donors<-nrow(D)/2;rr$residual_df<-nrow(D)-ncol(D)
  wt(rr,paste0(name,'_ROAST.tsv'));cat(name,': minimum FDR',min(rr$FDR,na.rm=TRUE),'\n');list(rr=rr,fit=fit)
 }
 fit<-run_roast(v,X,contrast,'primary',2026102607L)
 gt<-topTable(fit$fit,coef=ncol(X),n=Inf,sort.by='none',confint=.95)
 wt(data.frame(gene_id=rownames(gt),symbol=g$symbol[match(rownames(gt),g$gene_id)],gt,check.names=FALSE),'primary_gene_statistics.tsv')
 delta<-v$E[,sm$late==1,drop=FALSE]-v$E[,sm$late==0,drop=FALSE];colnames(delta)<-donors;wm(delta,'gene_paired_changes.tsv.gz')
 dd<-do.call(rbind,lapply(names(idx),function(gs)data.frame(gene_set=gs,donor=donors,paired_logCPM_change=colMeans(delta[idx[[gs]],,drop=FALSE]))));wt(dd,'program_donor_changes.tsv')
 flips<-as.matrix(expand.grid(rep(list(c(-1,1)),6)))
 ds<-do.call(rbind,lapply(names(idx),function(gs){z<-colMeans(delta[idx[[gs]],,drop=FALSE]);mu<-mean(z);se<-sd(z)/sqrt(6);ts<-as.vector(flips%*%z/6);data.frame(gene_set=gs,mean_change=mu,CI_L=mu-qt(.975,5)*se,CI_R=mu+qt(.975,5)*se,n_positive=sum(z>0),n_negative=sum(z<0),exact_signflip_p=mean(abs(ts)>=abs(mu)-1e-12))}));ds$signflip_FDR<-p.adjust(ds$exact_signflip_p,'BH',n=10);wt(ds,'program_donor_summary.tsv')
 loos<-list()
 for(i in seq_along(donors)){
  s<-sm[sm$donor!=donors[i],];j<-match(s$sample_column,samples)
  D<-X[j,colnames(X)!=donors[i],drop=FALSE];V<-normalize(C[keep,j,drop=FALSE],D)$v;ct<-c(rep(0,ncol(D)-1),1)
  z<-run_roast(V,D,ct,paste0('LOO_',sub('patient ','P',donors[i])),2026102640L+i);z$rr$omitted_donor<-donors[i];loos[[i]]<-z$rr
 }
 wt(do.call(rbind,loos),'all_LOO_ROAST.tsv')
 pc<-prcomp(t(v$E),center=TRUE,scale.=FALSE);wt(data.frame(sm,pc$x[,1:4],check.names=FALSE),'PCA_samples.tsv');wt(data.frame(PC=seq_along(pc$sdev),variance_fraction=pc$sdev^2/sum(pc$sdev^2)),'PCA_variance.tsv')
 safe_save(list(y=y,v=v,design=X,contrast=contrast,samples=sm,genes=g,members=members,indices=idx,fit=fit$fit),'analysis_objects.rds')
 status<-'COMPLETE';cat('All predefined longitudinal analyses complete.\n')
},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){error<<-conditionMessage(e);cat('ERROR',error,'\n')})
writeLines(warnings,file.path(OUT,'warnings.txt'));capture.output(sessionInfo(),file=file.path(OUT,'sessionInfo.txt'))
write_json(list(status=status,error=error,n_samples=12,n_donors=6,primary_model='paired late minus early',finished_utc=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(OUT,'RUN_STATUS.json'),auto_unbox=TRUE,pretty=TRUE)
if(status!='COMPLETE')stop(error)
}

analyse_radiology <- function() {
options(stringsAsFactors=FALSE,digits=16,width=180)
suppressPackageStartupMessages({library(edgeR);library(limma);library(jsonlite)})
ROOT <- 'data/airway/radiology'; OUT <- 'output/tables/airway/radiology'
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
rt<-function(p)read.delim(p,check.names=FALSE,comment.char='')
wt<-function(x,n)write.table(x,file.path(OUT,n),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
wm<-function(x,n){con<-gzfile(file.path(OUT,n),'wt');write.table(data.frame(gene_id=rownames(x),x,check.names=FALSE),con,sep='\t',row.names=FALSE,quote=TRUE);close(con);z<-rt(gzfile(file.path(OUT,n)));stopifnot(nrow(z)==nrow(x),ncol(z)==ncol(x)+1)}
safe_save<-function(x,n){con<-gzfile(file.path(OUT,n),'wb',compression=6);saveRDS(x,con);close(con);stopifnot(identical(readRDS(file.path(OUT,n)),x))}
RNGkind('Mersenne-Twister','Inversion','Rejection')
warnings<-character();status<-'FAILED';error<-NULL
tryCatch(withCallingHandlers({
 meta<-rt(file.path(ROOT,'input/early40_clinical_mapping.tsv'))
 all<-as.matrix(read.delim(gzfile(file.path(ROOT,'input/GSE261777_counts.txt.gz')),dec=',',check.names=FALSE,row.names=1))
 stopifnot(nrow(all)==34515,ncol(all)==55,all(is.finite(all)),all(all>=0),!anyDuplicated(rownames(all)))
 sm<-meta[!is.na(meta$PRA),];sm<-sm[order(sm$patient_id),];sm$sample_column<-as.character(sm$sample_column)
 samples<-sm$sample_column; stopifnot(nrow(sm)==37,sum(sm$PRA)==31,all(sm$visit=='study visit 1'),!anyDuplicated(samples))
 stopifnot(setequal(sm$patient_id[sm$PRA==0],c(4,6,7,19,21,36)))
 wt(sm,'clinical_sample_manifest.tsv');C<-all[,samples,drop=FALSE]
 X<-cbind(Intercept=1,PRA=sm$PRA);rownames(X)<-samples
 stopifnot(qr(X)$rank==2,nrow(X)-ncol(X)==35)
 ref<-rt(file.path('data/airway/assets/gencode_v44_symbols.tsv'));ref$stable<-sub('[.][0-9]+$','',ref$gene_id);ref$symbol[is.na(ref$symbol)]<-''
 owners<-split(ref$symbol,ref$stable);canon<-vapply(owners,function(s){u<-unique(s[nzchar(s)]);if(length(u)==1)u else ''},character(1))
 g<-data.frame(gene_id=rownames(C),symbol=unname(canon[match(rownames(C),names(canon))]));g$symbol[is.na(g$symbol)]<-'';g$mapped<-nzchar(g$symbol)
 y0<-DGEList(C);keep<-filterByExpr(y0,group=sm$PRA,min.count=10,min.total.count=15,large.n=10,min.prop=.7)
 g$primary_filter<-keep;g$mean_CPM<-rowMeans(cpm(y0));ord<-order(-g$mean_CPM,g$gene_id);ord<-ord[keep[ord]&g$mapped[ord]];reps<-ord[!duplicated(g$symbol[ord])]
 g$set_representative<-seq_len(nrow(g))%in%reps;wt(g,'gene_mapping_filter_details.tsv')
 panel<-rt(file.path('data/airway/assets/gene_set_panel.tsv'));members<-rt(file.path('data/airway/assets/gene_set_members.tsv'))
 hg<-rt(gzfile(file.path('data/airway/assets/hgnc_complete_set_20261017.txt.gz')));hg<-hg[hg$status=='Approved',]
 pr<-strsplit(ifelse(is.na(hg$prev_symbol),'',hg$prev_symbol),'|',fixed=TRUE);own<-split(rep(hg$symbol,lengths(pr)),unlist(pr,use.names=FALSE));prev<-vapply(own,function(s){u<-unique(s);if(length(u)==1)u else ''},character(1))
 members$symbol<-ifelse(members$gene_symbol%in%hg$symbol,members$gene_symbol,unname(prev[match(members$gene_symbol,names(prev))]));members$symbol[is.na(members$symbol)]<-''
 members$gene_id<-g$gene_id[reps][match(members$symbol,g$symbol[reps])];members$retained<-!is.na(members$gene_id);wt(members,'gene_set_member_details.tsv')
 idx<-setNames(lapply(panel$gene_set,function(gs)which(rownames(C)[keep]%in%members$gene_id[members$gene_set==gs&members$retained])),panel$gene_set)
 cov<-data.frame(gene_set=panel$gene_set,n_retained=lengths(idx),source_size=panel$set_size);cov$coverage<-cov$n_retained/cov$source_size;cov$eligible<-cov$n_retained>=10&cov$coverage>=.2;wt(cov,'program_coverage.tsv')
 normalize<-function(c,d){y<-calcNormFactors(DGEList(c),method='TMM');v<-voom(y,d,span=.5,normalize.method='none',plot=FALSE,save.plot=TRUE);dimnames(v$weights)<-dimnames(v$E);list(y=y,v=v)}
 run_model<-function(c,d,name,seed,save_gene=TRUE){
  stopifnot(ncol(c)==nrow(d),identical(colnames(c),rownames(d)),qr(d)$rank==ncol(d))
  nv<-normalize(c,d);v<-nv$v;ct<-c(rep(0,ncol(d)-1),1)
  set.seed(seed);rr<-mroast(v,index=idx[cov$eligible],design=d,contrast=ct,set.statistic='mean',nrot=19999L,robust=TRUE,trend=FALSE,approx.zscore=TRUE,legacy=FALSE,midp=FALSE,sort='none')
  rr<-data.frame(gene_set=rownames(rr),rr,check.names=FALSE);rr<-merge(cov,rr,by='gene_set',all.x=TRUE,sort=FALSE);rr<-rr[match(panel$gene_set,rr$gene_set),]
  rr$FDR<-p.adjust(rr$PValue,'BH',n=10);rr$FDR.Mixed<-p.adjust(rr$PValue.Mixed,'BH',n=10);rr$scenario<-name;rr$seed<-seed
  fit<-eBayes(lmFit(v,d),robust=TRUE,trend=FALSE)
  rr$mean_weighted_gene_logFC<-vapply(idx,function(i)mean(fit$coefficients[i,ncol(d)]),numeric(1));rr$n_samples<-nrow(d);rr$residual_df<-nrow(d)-ncol(d)
  wt(rr,paste0(name,'_ROAST.tsv'));wt(data.frame(sample=rownames(d),d,check.names=FALSE),paste0(name,'_design.tsv'))
  wt(data.frame(sample=colnames(c),nv$y$samples,check.names=FALSE),paste0(name,'_normalisation_factors.tsv'))
  if(save_gene){gt<-topTable(fit,coef=ncol(d),n=Inf,sort.by='none',confint=.95);wt(data.frame(gene_id=rownames(gt),symbol=g$symbol[match(rownames(gt),g$gene_id)],gt,check.names=FALSE),paste0(name,'_gene_statistics.tsv'))}
  cat(name,': minimum program FDR',min(rr$FDR,na.rm=TRUE),'\n');flush.console()
  list(y=nv$y,v=v,fit=fit,rr=rr,design=d)
 }
 primary<-run_model(C[keep,,drop=FALSE],X,'primary',2026102611L)
 wm(C,'clinical_submitted_counts.tsv.gz');wm(primary$v$E,'primary_normalised_logCPM.tsv.gz');wm(primary$v$weights,'primary_voom_weights.tsv.gz')
 score<-sapply(idx,function(ii)colMeans(primary$v$E[ii,,drop=FALSE]));rownames(score)<-samples
 wt(data.frame(sample=samples,patient_id=sm$patient_id,PRA=sm$PRA,score,check.names=FALSE),'program_sample_scores.tsv')
 frac_names<-colnames(rt(file.path(ROOT,'input/published_deconvolution_TableE3.tsv')))[2:17]
 pc_design<-function(s,name){
  f<-as.matrix(s[,frac_names]);stopifnot(all(is.finite(f)),all(f>=0),all(rowSums(f)>.99),all(rowSums(f)<1.01))
  f<-f/rowSums(f);h<-sqrt(f);pc<-prcomp(h,center=TRUE,scale.=FALSE)
  wt(data.frame(sample=s$sample_column,f,check.names=FALSE),paste0(name,'_fractions_closed.tsv'))
  wt(data.frame(sample=s$sample_column,pc$x,check.names=FALSE),paste0(name,'_composition_scores.tsv'))
  wt(data.frame(cell_type=rownames(pc$rotation),pc$rotation,check.names=FALSE),paste0(name,'_composition_loadings.tsv'))
  wt(data.frame(PC=seq_along(pc$sdev),variance_fraction=pc$sdev^2/sum(pc$sdev^2)),paste0(name,'_composition_variance.tsv'))
  D<-cbind(Intercept=1,PC1=pc$x[,1],PC2=pc$x[,2],PRA=s$PRA);rownames(D)<-as.character(s$sample_column)
  list(design=D,pc=pc,fractions=f)
 }
 pcs<-pc_design(sm,'composition37');composition<-run_model(C[keep,,drop=FALSE],pcs$design,'composition37',2026102612L)
 ss<-sm[sm$RMSE<.6,];stopifnot(nrow(ss)==34,sum(ss$PRA)==29,setequal(setdiff(sm$patient_id,ss$patient_id),c(9,23,36)))
 pcr<-pc_design(ss,'RMSE34_composition');restricted<-run_model(C[keep,ss$sample_column,drop=FALSE],pcr$design,'RMSE34_composition',2026102613L)
 loos<-list();controls<-sm$patient_id[sm$PRA==0]
 for(i in seq_along(controls)){
  s<-sm[sm$patient_id!=controls[i],];D<-X[as.character(s$sample_column),,drop=FALSE]
  z<-run_model(C[keep,as.character(s$sample_column),drop=FALSE],D,paste0('LOO_control_',controls[i]),2026102660L+i,FALSE);z$rr$omitted_control<-controls[i];loos[[i]]<-z$rr
 }
 wt(do.call(rbind,loos),'all_control_LOO_ROAST.tsv')
 clinical_genes<-rt(file.path(OUT,'primary_gene_statistics.tsv'));wt(clinical_genes[clinical_genes$symbol%in%c('SMURF1','HLA-DRA','NFKBIL1'),],'published_clinical_gene_details.tsv')
 lg<-rt('output/tables/airway/longitudinal/primary_gene_statistics.tsv');wt(lg[lg$symbol%in%c('PRSS3','HSPB8','HOPX','CAPN14','MMP12','MMP9'),],'published_longitudinal_gene_details.tsv')
 f<-primary$fit;wt(data.frame(gene_id=rownames(f),df_total=f$df.total,s2_post=f$s2.post,stdev_unscaled_PRA=f$stdev.unscaled[,2]),'primary_posterior_details.tsv')
 safe_save(list(primary=primary,composition=composition,restricted=restricted,samples=sm,genes=g,members=members,indices=idx,composition_full=pcs,composition_restricted=pcr,loo=loos),'clinical_analysis_objects.rds')
 status<-'COMPLETE';cat('All predefined clinical models complete.\n')
},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){error<<-conditionMessage(e);cat('ERROR',error,'\n')})
writeLines(warnings,file.path(OUT,'warnings.txt'));capture.output(sessionInfo(),file=file.path(OUT,'sessionInfo.txt'))
write_json(list(status=status,error=error,n_samples=37,n_with_PRA=31,n_without_PRA=6,primary_model='PRA_present_vs_absent',finished_utc=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(OUT,'RUN_STATUS.json'),auto_unbox=TRUE,pretty=TRUE)
if(status!='COMPLETE')stop(error)
}

analyse_pseudobulk <- function() {
options(stringsAsFactors=FALSE,digits=16,width=180)
suppressPackageStartupMessages({library(edgeR);library(limma);library(jsonlite)})
ROOT <- 'data/airway/pseudobulk'; OUT <- 'output/tables/airway/pseudobulk'
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
rt<-function(p)read.delim(p,check.names=FALSE,comment.char='')
wt<-function(x,n){dest<-file.path(OUT,n);tmp<-paste0(dest,'.part');write.table(x,tmp,sep='\t',row.names=FALSE,quote=TRUE,na='NA');back<-rt(tmp);stopifnot(nrow(back)==nrow(x),ncol(back)==ncol(x));stopifnot(file.rename(tmp,dest))}
wm<-function(x,n){dest<-file.path(OUT,n);tmp<-paste0(dest,'.part');con<-gzfile(tmp,'wt');write.table(data.frame(gene_id=rownames(x),x,check.names=FALSE),con,sep='\t',row.names=FALSE,quote=TRUE);close(con);back<-rt(gzfile(tmp));stopifnot(nrow(back)==nrow(x),ncol(back)==ncol(x)+1,all(is.finite(as.matrix(back[,-1]))));stopifnot(file.rename(tmp,dest))}
safe_save<-function(x,n){con<-gzfile(file.path(OUT,n),'wb',compression=6);saveRDS(x,con);close(con);stopifnot(identical(readRDS(file.path(OUT,n)),x))}
RNGkind('Mersenne-Twister','Inversion','Rejection')
warnings<-character();status<-'FAILED';error<-NULL
tryCatch(withCallingHandlers({
 sm<-rt(file.path(ROOT,'pseudobulk/donors.tsv'));rownames(sm)<-sm$donor
 celln<-rt(file.path(ROOT,'pseudobulk/donor_cell_counts.tsv'));rownames(celln)<-celln$donor
 covar<-rt(file.path(ROOT,'input/donor_covariates.tsv'));covar<-covar[match(sm$donor,covar$Study_ID),]
 covar$male<-as.integer(covar$Sex%in%c('M','Male'));covar$age_centered<-covar$Age-mean(covar$Age)
 stopifnot(nrow(sm)==11,sum(sm$case)==5,all(is.finite(covar$Age)),all(celln$primary>=50))
 wt(cbind(sm,covar[,c('published_ID','Age','male','covariate_linkage')]),'donor_manifest.tsv')
 get_counts<-function(name){z<-rt(gzfile(file.path(ROOT,'pseudobulk',paste0(name,'_counts.tsv.gz'))));m<-as.matrix(z[,-1]);rownames(m)<-z[[1]];storage.mode(m)<-'double';stopifnot(identical(colnames(m),sm$donor),!anyDuplicated(rownames(m)),all(is.finite(m)),all(m>=0),all(m==round(m)));m}
 C<-get_counts('primary')
 ref<-rt(file.path('data/airway/assets/gencode_v44_symbols.tsv'));ref$stable<-sub('[.][0-9]+$','',ref$gene_id);ref$symbol[is.na(ref$symbol)]<-''
 owners<-split(ref$symbol,ref$stable);canon<-vapply(owners,function(s){u<-unique(s[nzchar(s)]);if(length(u)==1)u else ''},character(1))
 genes<-data.frame(gene_id=rownames(C),symbol=unname(canon[match(rownames(C),names(canon))]));genes$symbol[is.na(genes$symbol)]<-''
 panel<-rt(file.path('data/airway/assets/gene_set_panel.tsv'));members0<-rt(file.path('data/airway/assets/gene_set_members.tsv'))
 hg<-rt(gzfile(file.path('data/airway/assets/hgnc_complete_set_20261017.txt.gz')));hg<-hg[hg$status=='Approved',]
 prevs<-strsplit(ifelse(is.na(hg$prev_symbol),'',hg$prev_symbol),'|',fixed=TRUE);owners<-split(rep(hg$symbol,lengths(prevs)),unlist(prevs,use.names=FALSE));prev<-vapply(owners,function(s){u<-unique(s);if(length(u)==1)u else ''},character(1))
 members0$symbol<-ifelse(members0$gene_symbol%in%hg$symbol,members0$gene_symbol,unname(prev[match(members0$gene_symbol,names(prev))]));members0$symbol[is.na(members0$symbol)]<-''
 mapping<-function(counts,s,name){
   y<-DGEList(counts);keep<-filterByExpr(y,group=s$case,min.count=10,min.total.count=15,large.n=10,min.prop=.7)
   g<-genes;g$retained<-keep;g$mean_CPM<-rowMeans(cpm(y));ord<-order(-g$mean_CPM,g$gene_id);ord<-ord[keep[ord]&nzchar(g$symbol[ord])];reps<-ord[!duplicated(g$symbol[ord])]
   g$set_representative<-seq_len(nrow(g))%in%reps
   m<-members0;m$gene_id<-g$gene_id[reps][match(m$symbol,g$symbol[reps])];m$retained<-!is.na(m$gene_id)
   idx<-setNames(lapply(panel$gene_set,function(gs)which(rownames(counts)[keep]%in%m$gene_id[m$gene_set==gs&m$retained])),panel$gene_set)
   cv<-data.frame(gene_set=panel$gene_set,n_retained=lengths(idx),source_size=panel$set_size);cv$coverage<-cv$n_retained/cv$source_size;cv$eligible<-cv$n_retained>=5
   wt(g,paste0(name,'_gene_mapping.tsv'));wt(m,paste0(name,'_member_mapping.tsv'));wt(cv,paste0(name,'_coverage.tsv'))
   list(keep=keep,genes=g,members=m,index=idx,coverage=cv)
 }
 mainmap<-mapping(C,sm,'primary')
 make_design<-function(s,adjusted=FALSE){
   if(adjusted){co<-covar[match(s$donor,covar$Study_ID),];d<-cbind(Intercept=1,age=co$age_centered,male=co$male,case=s$case)}
   else d<-cbind(Intercept=1,case=s$case)
   rownames(d)<-s$donor;d
 }
 run_model<-function(counts,s,name,seed,map,adjusted=FALSE,save_gene=TRUE){
   d<-make_design(s,adjusted);stopifnot(qr(d)$rank==ncol(d),identical(colnames(counts),s$donor))
   y<-calcNormFactors(DGEList(counts[map$keep,,drop=FALSE]),method='TMM')
   v<-voom(y,d,span=.5,normalize.method='none',plot=FALSE,save.plot=TRUE);dimnames(v$weights)<-dimnames(v$E)
   fit<-eBayes(lmFit(v,d),robust=TRUE,trend=FALSE)
   set.seed(seed);rr<-mroast(v,index=map$index[map$coverage$eligible],design=d,contrast=ncol(d),set.statistic='mean',nrot=19999L,robust=TRUE,trend=FALSE,approx.zscore=TRUE,legacy=FALSE,midp=FALSE,sort='none')
   rr<-data.frame(gene_set=rownames(rr),rr,check.names=FALSE);rr<-merge(map$coverage,rr,by='gene_set',all.x=TRUE,sort=FALSE);rr<-rr[match(panel$gene_set,rr$gene_set),]
   rr$PValue[!rr$eligible]<-1;rr$PValue.Mixed[!rr$eligible]<-1
   rr$FDR<-p.adjust(rr$PValue,'BH',n=10);rr$FDR.Mixed<-p.adjust(rr$PValue.Mixed,'BH',n=10)
   rr$mean_gene_logFC<-vapply(map$index,function(ii)if(length(ii))mean(fit$coefficients[ii,ncol(d)]) else NA_real_,numeric(1))
   rr$scenario<-name;rr$seed<-seed;rr$n_donors<-nrow(s);rr$n_case<-sum(s$case);rr$n_control<-sum(s$case==0);rr$residual_df<-nrow(d)-ncol(d);rr$design_condition_number<-kappa(d,exact=TRUE)
   wt(rr,paste0(name,'_ROAST.tsv'));wt(data.frame(donor=rownames(d),d,check.names=FALSE),paste0(name,'_design.tsv'))
   wt(data.frame(donor=s$donor,y$samples,check.names=FALSE),paste0(name,'_normalisation.tsv'))
   if(save_gene){gt<-topTable(fit,coef=ncol(d),n=Inf,sort.by='none',confint=.95);wt(data.frame(gene_id=rownames(gt),symbol=genes$symbol[match(rownames(gt),genes$gene_id)],gt,check.names=FALSE),paste0(name,'_gene_statistics.tsv'))}
   cat(name,'donors',nrow(s),'genes',nrow(v),'directional minimum q',min(rr$FDR),'\n');flush.console()
   list(y=y,v=v,fit=fit,results=rr,design=d,map=map,samples=s)
 }
 primary<-run_model(C,sm,'primary',2026102601L,mainmap)
 wm(primary$v$E,'primary_logCPM.tsv.gz');wm(primary$v$weights,'primary_weights.tsv.gz')
 wt(data.frame(gene_id=rownames(primary$fit),s2_post=primary$fit$s2.post,df_total=primary$fit$df.total,stdev_unscaled_case=primary$fit$stdev.unscaled[,2]),'primary_posterior_details.tsv')
 scores<-sapply(mainmap$index,function(ii)colMeans(primary$v$E[ii,,drop=FALSE]));wt(data.frame(donor=sm$donor,case=sm$case,scores,check.names=FALSE),'primary_program_scores.tsv')
 sens<-list()
 for(i in 1:2){nm<-c('no_mixed','strict_QC')[i];ss<-sm[celln[sm$donor,nm]>=50,,drop=FALSE];stopifnot(sum(ss$case)>=3,sum(ss$case==0)>=3);cc<-get_counts(nm)[,ss$donor,drop=FALSE];sens[[nm]]<-run_model(cc,ss,nm,2026102601L+i,mainmap)}
 sens[['age_sex']]<-run_model(C,sm,'age_sex',2026102604L,mainmap,TRUE)
 comp<-list();ctnames<-c('Basal','Differentiating','Ciliated','Secretory')
 for(i in seq_along(ctnames)){
   nm<-ctnames[i];ss<-sm[celln[sm$donor,nm]>=50,,drop=FALSE]
   if(sum(ss$case)<3||sum(ss$case==0)<3){wt(data.frame(scenario=nm,reason='fewer than 3 donors per group with 50 cells'),'skipped_compartment.tsv');next}
   cc<-get_counts(nm)[,ss$donor,drop=FALSE];mp<-mapping(cc,ss,nm);comp[[nm]]<-run_model(cc,ss,nm,2026102610L+i,mp)
 }
 cr<-do.call(rbind,lapply(comp,`[[`,'results'));cr$FDR_across_40<-p.adjust(cr$PValue,'BH',n=40);cr$FDR_Mixed_across_40<-p.adjust(cr$PValue.Mixed,'BH',n=40);wt(cr,'compartment_ROAST_FDR40.tsv')
 loo<-list()
 for(i in seq_len(nrow(sm))){ss<-sm[-i,,drop=FALSE];nm<-paste0('LOO_',sm$donor[i]);z<-run_model(C[,ss$donor,drop=FALSE],ss,nm,2026102620L+i,mainmap,FALSE,FALSE);z$results$omitted_donor<-sm$donor[i];loo[[i]]<-z$results}
 wt(do.call(rbind,loo),'all_donor_LOO_ROAST.tsv')
 safe_save(list(primary=primary,sensitivity=sens,compartments=comp,loo=loo),'donor_analysis_objects.rds')
 status<-'COMPLETE';cat('All planned donor analyses completed.\n')
},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){error<<-conditionMessage(e);cat('ERROR',error,'\n')})
writeLines(warnings,file.path(OUT,'warnings.txt'));capture.output(sessionInfo(),file=file.path(OUT,'sessionInfo.txt'))
write_json(list(status=status,error=error,analysis_unit='biological donor',contrast='case minus healthy',finished_utc=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(OUT,'RUN_STATUS.json'),auto_unbox=TRUE,pretty=TRUE)
if(status!='COMPLETE')stop(error)
}

analyse_blood <- function() {
options(stringsAsFactors=FALSE,digits=16,width=180)
suppressPackageStartupMessages({library(edgeR);library(limma);library(jsonlite)})
ROOT <- 'data/airway/blood'; OUT <- 'output/tables/airway/blood'
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
rt<-function(p)read.delim(p,check.names=FALSE,comment.char='')
wt<-function(x,n){dest<-file.path(OUT,n);tmp<-paste0(dest,'.part');write.table(x,tmp,sep='\t',row.names=FALSE,quote=TRUE,na='NA');b<-rt(tmp);stopifnot(nrow(b)==nrow(x),ncol(b)==ncol(x));stopifnot(file.rename(tmp,dest))}
wm<-function(x,n){dest<-file.path(OUT,n);tmp<-paste0(dest,'.part');con<-gzfile(tmp,'wt');write.table(data.frame(gene_id=rownames(x),x,check.names=FALSE),con,sep='\t',row.names=FALSE,quote=TRUE);close(con);con<-gzfile(tmp,'rt');b<-rt(con);close(con);stopifnot(nrow(b)==nrow(x),ncol(b)==ncol(x)+1,all(is.finite(as.matrix(b[,-1]))));stopifnot(file.rename(tmp,dest))}
safe_save<-function(x,n){dest<-file.path(OUT,n);tmp<-paste0(dest,'.part');con<-gzfile(tmp,'wb',compression=6);saveRDS(x,con);close(con);stopifnot(identical(readRDS(tmp),x));stopifnot(file.rename(tmp,dest))}
RNGkind('Mersenne-Twister','Inversion','Rejection');warnings<-character();status<-'FAILED';error<-NULL
tryCatch(withCallingHandlers({
 sm<-rt(file.path(ROOT,'input/primary_samples.tsv'));sm$donor<-sm$key;rownames(sm)<-sm$donor
 sm$age_centered<-sm$Age-mean(sm$Age);sm$log_time<-log2(sm$TPI)-mean(log2(sm$TPI))
 sm$burden<-sm$other_symptoms-mean(sm$other_symptoms)
 stopifnot(nrow(sm)==17,sum(sm$cough)==5,all(sm$TPI>=90),!anyDuplicated(sm$donor),sum(sm$male)==4)
 zz<-rt(gzfile(file.path(ROOT,'input/primary_counts.tsv.gz')));C<-as.matrix(zz[,-1]);rownames(C)<-zz[[1]];storage.mode(C)<-'double'
 stopifnot(identical(colnames(C),sm$donor),all(grepl('^ENSG[0-9]+$',rownames(C))),all(is.finite(C)),all(C>=0),all(C==round(C)))
 wt(sm,'sample_manifest.tsv')
 ref<-rt(file.path('data/airway/assets/gencode_v44_symbols.tsv'));ref$stable<-sub('[.][0-9]+$','',ref$gene_id);ref$symbol[is.na(ref$symbol)]<-''
 owners<-split(ref$symbol,ref$stable);canon<-vapply(owners,function(s){u<-unique(s[nzchar(s)]);if(length(u)==1)u else ''},character(1))
 genes<-data.frame(gene_id=rownames(C),symbol=unname(canon[match(rownames(C),names(canon))]));genes$symbol[is.na(genes$symbol)]<-''
 panel<-rt(file.path('data/airway/assets/gene_set_panel.tsv'));members<-rt(file.path('data/airway/assets/gene_set_members.tsv'))
 hg<-rt(gzfile(file.path('data/airway/assets/hgnc_complete_set_20261017.txt.gz')));hg<-hg[hg$status=='Approved',]
 prevs<-strsplit(ifelse(is.na(hg$prev_symbol),'',hg$prev_symbol),'|',fixed=TRUE);owners<-split(rep(hg$symbol,lengths(prevs)),unlist(prevs,use.names=FALSE));prev<-vapply(owners,function(s){u<-unique(s);if(length(u)==1)u else ''},character(1))
 members$symbol<-ifelse(members$gene_symbol%in%hg$symbol,members$gene_symbol,unname(prev[match(members$gene_symbol,names(prev))]));members$symbol[is.na(members$symbol)]<-''
 yy<-DGEList(C);keep<-filterByExpr(yy,group=sm$cough,min.count=10,min.total.count=15,large.n=10,min.prop=.7)
 genes$retained<-keep;genes$mean_CPM<-rowMeans(cpm(yy));ord<-order(-genes$mean_CPM,genes$gene_id);ord<-ord[keep[ord]&nzchar(genes$symbol[ord])];reps<-ord[!duplicated(genes$symbol[ord])]
 genes$set_representative<-seq_len(nrow(genes))%in%reps
 members$gene_id<-genes$gene_id[reps][match(members$symbol,genes$symbol[reps])];members$retained<-!is.na(members$gene_id)
 idx<-setNames(lapply(panel$gene_set,function(gs)which(rownames(C)[keep]%in%members$gene_id[members$gene_set==gs&members$retained])),panel$gene_set)
 cv<-data.frame(gene_set=panel$gene_set,n_retained=lengths(idx),source_size=panel$set_size);cv$coverage<-cv$n_retained/cv$source_size;cv$eligible<-cv$n_retained>=5
 wt(genes,'gene_mapping.tsv');wt(members,'member_mapping.tsv');wt(cv,'coverage.tsv')
 design<-function(s,mode){
   d<-switch(mode,primary=cbind(Intercept=1,age=s$age_centered,male=s$male,log2_TPI=s$log_time,cough=s$cough),
     unadjusted=cbind(Intercept=1,cough=s$cough),
     women_only=cbind(Intercept=1,age=s$age_centered,log2_TPI=s$log_time,cough=s$cough),
     other_symptoms=cbind(Intercept=1,age=s$age_centered,male=s$male,log2_TPI=s$log_time,other_symptoms=s$burden,cough=s$cough))
   rownames(d)<-s$donor;d
 }
 run<-function(s,name,seed,mode='primary',save_gene=TRUE){
   d<-design(s,mode);stopifnot(qr(d)$rank==ncol(d),nrow(d)>ncol(d),identical(rownames(d),s$donor))
   y<-calcNormFactors(DGEList(C[keep,s$donor,drop=FALSE]),method='TMM')
   v<-voom(y,d,span=.5,normalize.method='none',plot=FALSE,save.plot=TRUE);dimnames(v$weights)<-dimnames(v$E)
   fit<-eBayes(lmFit(v,d),robust=TRUE,trend=FALSE)
   set.seed(seed);rr<-mroast(v,index=idx[cv$eligible],design=d,contrast=ncol(d),set.statistic='mean',nrot=19999L,robust=TRUE,trend=FALSE,approx.zscore=TRUE,legacy=FALSE,midp=FALSE,sort='none')
   rr<-data.frame(gene_set=rownames(rr),rr,check.names=FALSE);rr<-merge(cv,rr,by='gene_set',all.x=TRUE,sort=FALSE);rr<-rr[match(panel$gene_set,rr$gene_set),]
   rr$PValue[!rr$eligible]<-1;rr$PValue.Mixed[!rr$eligible]<-1;rr$FDR<-p.adjust(rr$PValue,'BH',n=10);rr$FDR.Mixed<-p.adjust(rr$PValue.Mixed,'BH',n=10)
   rr$mean_gene_logFC<-vapply(idx,function(ii)if(length(ii))mean(fit$coefficients[ii,ncol(d)]) else NA_real_,numeric(1))
   rr$scenario<-name;rr$seed<-seed;rr$n_donors<-nrow(s);rr$n_cough<-sum(s$cough);rr$n_not_listed<-sum(s$cough==0);rr$residual_df<-nrow(d)-ncol(d);rr$design_condition_number<-kappa(d,exact=TRUE)
   wt(rr,paste0(name,'ROAST.tsv'));wt(data.frame(donor=rownames(d),d,check.names=FALSE),paste0(name,'design.tsv'));wt(data.frame(donor=s$donor,y$samples,check.names=FALSE),paste0(name,'normalisation.tsv'))
   gt<-topTable(fit,coef=ncol(d),n=Inf,sort.by='none',confint=.95)
   if(save_gene)wt(data.frame(gene_id=rownames(gt),symbol=genes$symbol[match(rownames(gt),genes$gene_id)],gt,check.names=FALSE),paste0(name,'gene_statistics.tsv'))
   cat(name,'n',nrow(s),'df',nrow(d)-ncol(d),'genes',nrow(v),'minimum directional q',min(rr$FDR),'\n');flush.console()
   list(y=y,v=v,fit=fit,results=rr,design=d,samples=s,genes=gt)
 }
 primary<-run(sm,'primary_',2026102701L)
 wm(primary$v$E,'primary_logCPM.tsv.gz');wm(primary$v$weights,'primary_weights.tsv.gz')
 wt(data.frame(gene_id=rownames(primary$fit),s2_post=primary$fit$s2.post,df_total=primary$fit$df.total,stdev_unscaled_cough=primary$fit$stdev.unscaled[,ncol(primary$design)]),'primary_posterior_details.tsv')
 sc<-sapply(idx,function(ii)colMeans(primary$v$E[ii,,drop=FALSE]));wt(data.frame(donor=sm$donor,cough=sm$cough,sc,check.names=FALSE),'program_scores.tsv')
 sens<-list(unadjusted=run(sm,'unadjusted_',2026102702L,'unadjusted'),
             women_only=run(sm[sm$male==0,],'women_only_',2026102703L,'women_only'),
             other_symptoms=run(sm,'other_symptoms_',2026102704L,'other_symptoms'))
 loo<-list();lfc<-matrix(NA_real_,sum(keep),nrow(sm),dimnames=list(rownames(C)[keep],sm$donor));lq<-lfc
 for(i in seq_len(nrow(sm))){
   ss<-sm[-i,];dd<-design(ss,'primary')
   if(qr(dd)$rank<ncol(dd)){wt(data.frame(omitted=sm$donor[i],reason='rank deficient'),paste0('LOO_',sm$donor[i],'_SKIPPED.tsv'));next}
   z<-run(ss,paste0('LOO_',sm$donor[i],'_'),2026102720L+i,'primary',FALSE)
   z$results$omitted_donor<-sm$donor[i];loo[[sm$donor[i]]]<-z$results
   lfc[,i]<-z$genes$logFC;lq[,i]<-z$genes$adj.P.Val
 }
 wt(do.call(rbind,loo),'all_LOO_ROAST.tsv');wm(lfc,'LOO_gene_logFC.tsv.gz');wm(lq,'LOO_gene_FDR.tsv.gz')
 safe_save(list(primary=primary,sensitivity=sens,loo=loo,index=idx,coverage=cv,gene_mapping=genes,members=members),'cough_annotation_models.rds')
 status<-'COMPLETE';cat('All planned cough-annotation models completed.\n')
},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){error<<-conditionMessage(e);cat('ERROR',error,'\n')})
writeLines(warnings,file.path(OUT,'warnings.txt'));capture.output(sessionInfo(),file=file.path(OUT,'sessionInfo.txt'))
write_json(list(status=status,error=error,endpoint='cough_recorded_vs_not_listed',unit='one participant',finished_utc=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(OUT,'RUN_STATUS.json'),auto_unbox=TRUE,pretty=TRUE)
if(status!='COMPLETE')stop(error)
}

exact_scores <- function(scores, case, tolerance) {
    scores <- as.matrix(scores); n <- nrow(scores); n1 <- sum(case); n0 <- n - n1
    stopifnot(n1 > 1, n0 > 1, all(is.finite(scores)))
    a <- scores[case, , drop=FALSE]; b <- scores[!case, , drop=FALSE]
    difference <- colMeans(a) - colMeans(b)
    va <- apply(a, 2, var)/n1; vb <- apply(b, 2, var)/n0
    se <- sqrt(va + vb); df <- (va + vb)^2/(va^2/(n1-1) + vb^2/(n0-1))
    observed <- difference/se
    k <- min(n0, n1); indices <- combn(n, k)
    total <- colSums(scores); square <- colSums(scores^2)
    extreme <- integer(ncol(scores)); count <- ncol(indices)
    for (start in seq.int(1L, count, by=16384L)) {
        ind <- indices[, start:min(count, start+16383L), drop=FALSE]
        for (j in seq_len(ncol(scores))) {
            z <- matrix(scores[ind, j], nrow=k)
            s <- colSums(z); q <- colSums(z^2)
            other <- total[j]-s; otherq <- square[j]-q
            v <- pmax(0, (q-s^2/k)/(k-1))
            vo <- pmax(0, (otherq-other^2/(n-k))/(n-k-1))
            den <- sqrt(v/k+vo/(n-k))
            tt <- ifelse(den > 0, (s/k-other/(n-k))/den, 0)
            extreme[j] <- extreme[j]+sum(abs(tt) >= abs(observed[j])-tolerance)
        }
    }
    p <- extreme/count
    data.frame(gene_set=colnames(scores), difference=difference, SE=se,
               Welch_df=df, CI_L=difference-qt(.975,df)*se,
               CI_R=difference+qt(.975,df)*se, Welch_t=observed,
               exact_allocations=count, extreme_allocations=extreme,
               exact_P=p, exact_FDR=p.adjust(p,'BH'), row.names=NULL)
}

summarise_airway <- function() {
    base <- 'output/tables/airway'
    read <- function(path) read.delim(path, check.names=FALSE)
    write <- function(x, path) write.table(x, path, sep='\t', row.names=FALSE, quote=TRUE, na='NA')
    for (kind in c('radiology', 'pseudobulk', 'blood')) {
        name <- switch(kind, radiology='program_sample_scores.tsv',
                       pseudobulk='primary_program_scores.tsv', blood='program_scores.tsv')
        x <- read(file.path(base,kind,name))
        indicator <- switch(kind,radiology='PRA',pseudobulk='case',blood='cough')
        panel <- read('data/airway/assets/gene_set_panel.tsv')$gene_set
        z <- exact_scores(x[,panel,drop=FALSE], x[[indicator]]==1,
                          if(kind=='radiology') 1e-10 else 1e-12)
        stopifnot(unique(z$exact_allocations)==switch(kind,radiology=2324784,pseudobulk=462,blood=6188))
        write(z,file.path(base,kind,'score_exact_permutation.tsv'))
    }
    cells <- read(gzfile('data/airway/pseudobulk/pseudobulk/cell_QC_and_annotations.tsv.gz'))
    manifest <- read(file.path(base,'pseudobulk/donor_manifest.tsv'))
    tab <- table(cells$donor,cells$cell_type)
    proportions <- prop.table(tab,1)[manifest$donor,,drop=FALSE]
    z <- exact_scores(proportions,manifest$case==1,1e-12)
    names(z)[names(z)=='gene_set'] <- 'cell_type'
    write(z,file.path(base,'pseudobulk/cell_composition_exact_permutation.tsv'))
    rows <- lapply(c('culture','longitudinal','radiology','pseudobulk','blood'), function(kind) {
        z <- read(file.path(base,kind,'primary_ROAST.tsv'))
        effect <- if('mean_weighted_gene_logFC'%in%names(z)) z$mean_weighted_gene_logFC else z$mean_gene_logFC
        data.frame(context=kind,gene_set=z$gene_set,directional_P=z$PValue,
                   directional_q=z$FDR,mixed_P=z$PValue.Mixed,mixed_q=z$FDR.Mixed,
                   mean_gene_logFC=effect)
    })
    write(do.call(rbind,rows),file.path(base,'primary_programmes.tsv'))
}

analyses <- list(culture=analyse_culture, longitudinal=analyse_longitudinal,
                radiology=analyse_radiology, pseudobulk=analyse_pseudobulk, blood=analyse_blood)
selected <- commandArgs(trailingOnly=TRUE)
if (!length(selected)) selected <- names(analyses)
stopifnot(all(selected %in% c(names(analyses),'summaries')))
for (name in intersect(selected,names(analyses))) {
    message('Analysing ',name)
    analyses[[name]]()
}
if (all(names(analyses) %in% selected) || 'summaries' %in% selected) summarise_airway()
