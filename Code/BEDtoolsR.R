#### BEDtools R functions

bedTools.2in<-function(functionstring="/Users/Gabriel/Desktop/bin/bedtools2/bin/intersectBed",bed1,bed2,opt.string=""){
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen=99) # not to use scientific notation when writing out

  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  file.info(out)$size->size
  if(size>0){	
  	res=read.table(out,header=F)
  	unlink(a.file);unlink(b.file);unlink(out)
  	return(res)
  }	
}

bedTools.2merge<-function(functionstring="/Users/Gabriel/Desktop/bin/bedtools2/bin/mergeBed",bed1,opt.string=""){
  #create temp files
  a.file=tempfile()
  out   =tempfile()
  options(scipen =99)
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)

  command=paste(functionstring,"-i",a.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=read.table(out,header=F)
  unlink(a.file);unlink(out)
  return(res)
}
bedTools.2sort<-function(functionstring="/Users/Gabriel/Desktop/bin/bedtools2/bin/sortBed",bed1,opt.string=""){
  #create temp files
  a.file=tempfile()
  out   =tempfile()
  options(scipen =99)
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  command=paste(functionstring,"-i",a.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=read.table(out,header=F)
  unlink(a.file);unlink(out)
  return(res)
}

bedTools.2closest<-function(functionstring="/Users/Gabriel/Desktop/bin/bedtools2/bin/closestBed",bed1,bed2,opt.string=""){
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen=99) # not to use scientific notation when writing out

  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  file.info(out)$size->size
  if(size>0){	
  	res=read.table(out,header=F)
  	unlink(a.file);unlink(b.file);unlink(out)
  	return(res)
  }	
}


bedTools.2sub<-function(functionstring="/Users/Gabriel/Desktop/bin/bedtools2/bin/subtractBed",bed1,bed2,opt.string=""){
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen=99) # not to use scientific notation when writing out

  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  file.info(out)$size->size
  if(size>0){	
  	res=read.table(out,header=F)
  	unlink(a.file);unlink(b.file);unlink(out)
  	return(res)
  }	
}