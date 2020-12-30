
d<- dir("~/Desktop/BAM_to_csv", full.names=TRUE)
newcall = NULL

for(x in 1:length(d)) {
		#count.fields(d[x])
		dat = read.table(d[x], header=T, comment.char="", stringsAsFactors=F, sep="\t")
		newcall = rbind(newcall, dat, deparse.level=0)
		print(paste(d[x],": \n", dim(dat)))	
}

d<-read.csv("/Users/rr11/Documents/SFHS_Project/De_novo_validation/BAM_to_csv/603_159848_output.csv", header=T, colClasses="character")

f=newcall
 
  Child_Tot<-(f$child_A + f$child_T + f$child_G + f$child_C + f$child_N + f$child_indel_ref + f$child_indel_alt + f$child_indel_other)  
  Father_Tot<- (f$father_A + f$father_T + f$father_G + f$father_C + f$father_N +f$father_indel_ref + f$father_indel_alt + f$father_indel_other) 
  Mother_Tot<- (f$mother_A + f$mother_T + f$mother_G + f$mother_C + f$mother_N +f$mother_indel_ref + f$mother_indel_alt + f$mother_indel_other)
 
  f<-cbind(f,Child_Tot,Father_Tot,Mother_Tot)
  
  
Child_ALT=vector() 
for (j in 1:length(f[,1])){
  
  		if (nchar(as.character(f$ALT[j]))!=1 |nchar(as.character(f$REF[j]))!=1) {				# Test whether one of the alleles has more than one base
			Child_ALT[j] <- as.character(f$child_indel_alt[j]+ f$child_indel_other[j])					
		} 
		
		if (nchar(as.character(f$ALT[j]))==1 & nchar(as.character(f$REF[j]))==1) {
		
			if (f$ALT[j]=="T"){
  				Child_ALT[j]<- as.character(f$child_T[j])
  			}
  			if (f$ALT[j]=="C"){
  				Child_ALT[j]<- as.character(f$child_C[j])
  			}
  			if (f$ALT[j]=="G"){
  				Child_ALT[j]<- as.character(f$child_G[j])
  			}
  			if (f$ALT[j]=="A"){
  				Child_ALT[j]<- as.character(f$child_A[j])
  			}
  		} 	
	}
 Paternal_ALT=vector() 
for (j in 1:length(f[,1])){
  
  		if (nchar(as.character(f$ALT[j]))!=1 | nchar(as.character(f$REF[j]))!=1) {				# Test whether one of the alleles has more than one base
			Paternal_ALT[j] <- as.character(f$father_indel_alt[j] + f$father_indel_other[j])					
		} 
		
		if (nchar(as.character(f$ALT[j]))==1 & nchar(as.character(f$REF[j]))==1) {
			if (f$ALT[j]=="T"){
  				Paternal_ALT[j]<- as.character(f$father_T[j])
  			}
  			if (f$ALT[j]=="C"){
  				Paternal_ALT[j]<- as.character(f$father_C[j])
  			}
  			if (f$ALT[j]=="G"){
  				Paternal_ALT[j]<- as.character(f$father_G[j])
  			}
  			if (f$ALT[j]=="A"){
  				Paternal_ALT[j]<- as.character(f$father_A[j])
  			}
		}	
  }
  Maternal_ALT=vector() 
for (j in 1:length(f[,1])){
  
  		if (nchar(as.character(f$ALT[j]))!=1 |nchar(as.character(f$REF[j]))!=1) {				# Test whether one of the alleles has more than one base
			Maternal_ALT[j] <- as.character(f$mother_indel_alt[j] + f$mother_indel_other[j])					
		} 
		
		if (nchar(as.character(f$ALT[j]))==1 & nchar(as.character(f$REF[j]))==1) {
			if (f$ALT[j]=="T"){
  				Maternal_ALT[j]<- as.character(f$mother_T[j])
  			}
  			if (f$ALT[j]=="C"){
  				Maternal_ALT[j]<- as.character(f$mother_C[j])
  			}
  			if (f$ALT[j]=="G"){
  				Maternal_ALT[j]<- as.character(f$mother_G[j])
  			}
  			if (f$ALT[j]=="A"){
  				Maternal_ALT[j]<- as.character(f$mother_A[j])
  			}
		} 
  	} 
  
  
f<-cbind(f, Child_ALT, Paternal_ALT, Maternal_ALT)
  
write.table(f, file="/Users/rr11/Documents/SFHS_Project/De_novo_validation/new_BAM_to_txt/Miseq_outpot_per_trio/Miseq_output_244_73430.txt", sep="\t", row.names=F, col.names=T, quote=F)

  
# example validation code for resequencing of candidate DNMs

num.candidates<-length(candidate.dnms[,1])

 #asigning the empty values to the models
	
#mut.m=as.numeric(mut.m)
#mut.d=as.numeric(mut.d)
#mut.c=as.numeric(mut.c)

	if(length(models[which(LL.models==max(LL.models))])==1) { 
	
		
		} else {best.model[i]<-NA}
	
	best.model.filt<-best.model
	
	best.model.filt[which(LRT<5)]<-NA
	


  
  
  