### ################################################################################################
setwd("~/Documents/BU/2017SPRING/Characteristic_Respondents/demo")

### load some necessary pacakges
library(sem)                            ## read in the SEM package
library(psych)                          ## loads package for cronbach alpha


### ################################################################################################
### read in the data...the master data frame (MDF)
### no need for cleaning or reformatting
###
MDF <- read.csv("sem_congeneric_demo_data.csv",stringsAsFactors=FALSE)


### ################################################################################################
### automate the SEM model building process
###
### this was the first pass thru all of the congeneric models
###

### ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
### create lists for the indices for items for each factor
inds.list <- list(c(1:5),
                  c(1:6))

# sink(file="tmp_output_170130v1.txt")
### ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
### use lists to build & run the congeneric models
for(i in 1:length(inds.list)) {
   cat(paste(c("\n\n\n••••••••••\n••• i = ",i,"\n·········\n\n"),collapse=""))
   ### ················································
   ### build the SEM models
   vars.4.model <- names(MDF)[inds.list[[i]]]
   tmp.text <- vars.4.model
   tmp.text <- paste("   PIPS_F1 --> ",tmp.text,sep="")
   tmp.text <- paste(tmp.text,", lambda",sep="")
   tmp.text <- paste(tmp.text,1:length(vars.4.model),sep="")
   tmp.text <- paste(tmp.text,"1,  NA\n",sep="")
   tmp.text <- paste(tmp.text,collapse="")
   tmp.text <- paste("   PIPS_F1     <-> PIPS_F1,               NA,        1\n",tmp.text,collapse="")
   cat(tmp.text)
   model.sem <- specifyModel(text=tmp.text,quiet=TRUE)

   eval(parse(text=paste(c("model.sem.",i," <- model.sem"),collapse="")))

   ### ················································
   ### run the SEM models
   use.vars <- vars.4.model
   tmp.data <- na.omit(MDF[,use.vars])
   NN <- dim(tmp.data)[1]
   S.cov.mat <- cov(tmp.data)
   sem.out <- sem(model.sem,S.cov.mat,NN,maxiter=10^5)
   eval(parse(text=paste(c("sem.out.",i," <- sem.out"),collapse="")))

   ### ················································
   ### generate the output for the SEM models
   print(summary(sem.out,fit.indices = c("RMSEA", "NNFI", "CFI", "SRMR")))
   print(stdCoef(sem.out))
   sem.modin <- modIndices(sem.out)
   n.for.A <- length(which(na.omit(as.vector(sem.modin$mod.A)) > 9))
   n.for.P <- length(which(na.omit(as.vector(sem.modin$mod.P)) > 9))/2
   print(sem.modin,n.largest=max(n.for.A,n.for.P))
   }
# sink()

rm(tmp.text,model.sem,sem.out,vars.4.model)
rm(i,use.vars,tmp.data,NN,S.cov.mat,sem.modin,n.for.A,n.for.P)


add.list <- list(c("Q02<->Q05","err0205",NA),
                 c("Q02<->Q05","err0205",NA))



# sink(file="tmp_output_170130v2.txt")
### ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
### use lists to build & run the congeneric models
for(i in 1:length(inds.list)) {
   cat(paste(c("\n\n\n••••••••••\n••• i = ",i,"\n·········\n\n"),collapse=""))
   ### ················································
   ### build the SEM models
   vars.4.model <- names(MDF)[inds.list[[i]]]
   tmp.text <- vars.4.model
   tmp.text <- paste("   PIPS_F1 --> ",tmp.text,sep="")
   tmp.text <- paste(tmp.text,", lambda",sep="")
   tmp.text <- paste(tmp.text,1:length(vars.4.model),sep="")
   tmp.text <- paste(tmp.text,"1,  NA\n",sep="")
   tmp.text <- paste(tmp.text,collapse="")
   tmp.text <- paste("   PIPS_F1     <-> PIPS_F1,               NA,        1\n",tmp.text,collapse="")
   # cat(tmp.text)
   model.sem <- specifyModel(text=tmp.text,quiet=TRUE)

if(typeof(add.list[[i]])=="character") {
   new.row <- add.list[[i]]
   model.sem <- rbind(model.sem,new.row)
   }

   eval(parse(text=paste(c("model.sem.",i," <- model.sem"),collapse="")))

   ### ················································
   ### run the SEM models
   use.vars <- vars.4.model
   tmp.data <- na.omit(MDF[,use.vars])
   NN <- dim(tmp.data)[1]
   S.cov.mat <- cov(tmp.data)
   sem.out <- sem(model.sem,S.cov.mat,NN,maxiter=10^5)
   eval(parse(text=paste(c("sem.out.",i," <- sem.out"),collapse="")))

   ### ················································
   ### generate the output for the SEM models
   print(summary(sem.out,fit.indices = c("RMSEA", "NNFI", "CFI", "SRMR")))
   print(stdCoef(sem.out))
   sem.modin <- modIndices(sem.out)
   n.for.A <- length(which(na.omit(as.vector(sem.modin$mod.A)) > 9))
   n.for.P <- length(which(na.omit(as.vector(sem.modin$mod.P)) > 9))/2
   print(sem.modin,n.largest=max(n.for.A,n.for.P))
   }
# sink()

rm(tmp.text,model.sem,sem.out,vars.4.model)
rm(i,use.vars,tmp.data,NN,S.cov.mat,sem.modin,n.for.A,n.for.P)









# factor.scores <- fscores(sem.out,tmp.data)
# avg.scores <- apply(tmp.data,1,FUN="mean",na.rm=T)
# plot(avg.scores,factor.scores[,1],pch=16,cex=.4)
# 1-cor(avg.scores,factor.scores[,1])^2


