nurigCalc <- function(reference, query){
  library(dplyr)
  library(tidyr)
  
  
 
  
  ref = reference$datapath
  lista = query
  
  
 
  
  # if (ncol(lista) == 2)
  # {
  #   colnames(lista) = c("file", "names")
  # } else{
  #   colnames(lista) = c("file")
  #   lista$names <- lista$file
  # }
  lista = rename(lista, file = datapath, names = name)
  
  
  j = 1
  for (i in lista$file)
  {
    #system(paste("./bin/mummer/nucmer ", ref, i, collapse = ""))
    #system("./bin/mummer/delta-filter -q out.delta > out.f.delta ")
    #system("./bin/mummer/show-coords -l -L 200 -q -B out.f.delta > out.coords")
    system(paste("nucmer ", ref, i, collapse = ""))
    system("delta-filter -q out.delta > out.f.delta ")
    system("show-coords -l -L 200 -q -B out.f.delta > out.coords")
    
    tmpTable = read.table("out.coords", sep = "\t")
    colnames(tmpTable) = c(
      "qID",
      "date",
      "lenq",
      "alitype",
      "referencefile",
      "rID",
      "startquery",
      "endquery",
      "startreference",
      "endreference",
      "peridentity",
      "persimilarity",
      "lenAlig",
      "c16",
      "c17",
      "c18",
      "c19",
      "strand",
      "lenref",
      "c21",
      "c22"
    )
    tmpTable$slot = j
    tmpTable$file = i
    
    if (j > 1)
    {
      finalTable = bind_rows(finalTable, tmpTable)
    } else{
      finalTable = tmpTable
    }
    j = j + 1
  }
  
  finalTable = full_join(finalTable, lista, by = "file")
  
  finalTable = finalTable %>% group_by(rID) %>% summarise(ln = first(lenref)) %>% mutate(prueba = cumsum(ln)) %>% mutate(position = prueba -
                                                                                                                           ln) %>% select(rID, position) %>% full_join(finalTable, .)
  finalTable = finalTable %>% group_by(slot) %>% mutate(orden = sum(lenAlig *
                                                                      peridentity / 100)) %>% arrange(desc(orden))
  
  return(finalTable)
}

nurigRender <- function(finalTable,reference,pal, formatImage, height, width, annotation, binPath){
  library(dplyr)
  library(tidyr)
  library(XML)
  
  if(missing(pal))
  {
    palette = 'rainbow'
  }
  if(missing(formatImage))
  {
    formatImage = 'SVG'
  }
  if(missing(height))
  {
    height = 800
  }
  
  if(missing(width))
  {
    width = 800
  }
  
  outFile = "out.xml"
  
  if(missing(annotation) | is.null(annotation))
  {
    annot = NULL
  }else{
    annot <- annotation
  }
  palette = pal
  
  tmp = finalTable %>%  split(.$slot)
  size = finalTable %>% group_by(rID) %>% summarise(tl = first(lenref)) %>% summarise(S = sum(tl))
  xml = xmlTree()
  xml$addTag(
    "cgview",
    attrs = c(
      backboneRadius = "160",
      sequenceLength = as.integer(size[1]),
      title = reference$name,
      labelPlacementQuality = "best"
    ),
    close = FALSE
  )
  
  
  
  j = 1
  
  if (palette == 'rainbow')
  {
    colores = rainbow(max(finalTable$slot))
  } else if (palette == 'topo') {
    colores = topo.colors(max(finalTable$slot))
  } else if (palette == 'terrain') {
    colores = terrain.colors(max(finalTable$slot))
  }
  
  for (tab in 1:length(tmp))
  {
    xml$addTag("featureSlot",
               attrs = c(strand = "direct"),
               close = FALSE)
    clr = paste("rgb(", paste(col2rgb(colores[j]), collapse = ","), ")", collapse = "")
    j = j + 1
    xml$addTag("feature",
               attrs = c(color = clr, decoration = "arc"),
               close = FALSE)
    for (r in 1:nrow(tmp[[tab]]))
    {
      #xml$addTag("feature", close = FALSE)
      xml$addTag(
        "featureRange",
        attrs = c(
          start = tmp[[tab]]$startreference[r] + tmp[[tab]]$position[r],
          stop = tmp[[tab]]$endreference[r] + tmp[[tab]]$position[r],
          opacity = (tmp[[tab]]$peridentity[r] / 100)^2
        ),
        close = FALSE
      )
      xml$closeTag()
      
    }
    xml$closeTag()
    xml$closeTag()
    
    
  }
  xml$closeTag()
  
  ####Draw Legend
  
  xml$addTag("legend",
             attrs = c(position = "upper-right"),
             close = FALSE)
  j = 1
  
  tmp2 = finalTable %>% group_by(file, names, orden) %>% summarise(slot = first(slot)) %>% arrange(desc(orden))
  
  
  
  for (i in tmp2$names)
  {
    clr = paste("rgb(", paste(col2rgb(colores[j]), collapse = ","), ")", collapse = "")
    xml$addTag(
      "legendItem",
      attrs = c(
        text = i,
        drawSwatch = "true",
        swatchColor = clr
      ),
      close = TRUE
    )
    j = j + 1
  }
  xml$closeTag()
  
  ##### Draw contigs
  
  contigs = finalTable %>% group_by(rID) %>% summarise(position = first(position) +
                                                         1)
  lblBoolean = (!is.null(annot) && nrow(annot)>0 && ncol(annot)< 9)
  
  if (nrow(contigs) > 1)
  {
    xml$addTag("featureSlot",
               attrs = c(strand = "direct"),
               close = FALSE)
    xml$addTag("feature",
               attrs = c(decoration = "clockwise-arrow"),
               close = FALSE)
    cl = rep(c("black", "grey"), times = (nrow(contigs) / 2) + 1)
    for (i in 1:(nrow(contigs) - 1))
    {
      xml$addTag(
        "featureRange",
        attrs = c(
          start = contigs$position[i],
          stop = contigs$position[i + 1],
          mouseover = contigs$rID[i],
          label = contigs$rID[i],
          color = cl[i],
          showLabel = !lblBoolean
        ),
        close = TRUE
      )
    }
    xml$addTag(
      "featureRange",
      attrs = c(
        start = contigs$position[nrow(contigs)]
        ,
        stop = as.integer(size[1]),
        mouseover = contigs$rID[nrow(contigs)],
        label = contigs$rID[nrow(contigs)],
        color = cl[i + 1],
        showLabel = !lblBoolean
      ),
      close = TRUE
    )
    xml$closeTag()
    xml$closeTag()
  }
  
  ### Draw Annotation
  
  if(!is.null(annot) && nrow(annot)>0 && ncol(annot)< 9)
  {
    
    annot = annot %>% rename(rID = seqname)
    
    
    #colnames(annot) = c("rID","start","end","strand","label","color")
    finalAnnot <- finalTable %>% group_by(rID) %>% summarise(position = first(position)) %>% inner_join(., annot) %>%  mutate(start = start +
                                                                                                                               position, end = end + position)
    
    xml$addTag("featureSlot",
               attrs = c(strand = "direct"),
               close = FALSE)
    xml$addTag("feature",
               attrs = c(decoration = "clockwise-arrow"),
               close = FALSE)
    cl = rep(c("red", "blue"), times = (nrow(finalAnnot) / 2) + 1)
    tmp = finalAnnot %>% filter(strand == "+")
    for (i in 1:(nrow(tmp)))
    {
      xml$addTag(
        "featureRange",
        attrs = c(
          start = tmp$start[i],
          stop = tmp$end[i],
          color = tmp$color[i],
          label = tmp$label[i],
          showLabel = as.character(!(is.na(tmp$label[i])|| tmp$label[i]==""))
        ),
        close = TRUE
      )
    }
    xml$closeTag()
    xml$closeTag()
    xml$addTag("featureSlot",
               attrs = c(strand = "direct"),
               close = FALSE)
    xml$addTag("feature",
               attrs = c(decoration = "counterclockwise-arrow"),
               close = FALSE)
    tmp = finalAnnot %>% filter(strand == "-")
    for (i in 1:(nrow(tmp)))
    {
      xml$addTag(
        "featureRange",
        attrs = c(
          start = tmp$start[i],
          stop = tmp$end[i],
          color = tmp$color[i],
          label = tmp$label[i],
          showLabel = as.character(!(is.na(tmp$label[i])|| tmp$label[i]==""))
        ),
        close = TRUE
      )
    }
    xml$closeTag()
    xml$closeTag()
  }
  #binPath = gsub("[[:space:]]","",binPath)
  saveXML(xml$value(), file = "out.xml")
  system(paste("java -jar ",binPath,"/bin/cgview/cgview.jar -i out.xml -o fig.svg -f svg -W ",width," -H ",height,collapse ="", sep =""))
  #write.table(finalTable,file = "FinalTable.tsv",sep = '\t', row.names = FALSE)
}

gffImporter <- function(path)
{
  gff = read.delim(as.character(path), sep = "\t",comment.char = "#", header = FALSE, stringsAsFactors = FALSE)
  colnames(gff) = c("seqname","source","feature","start","end","score","strand","frame","attribute")
  gff$Key = seq(1:nrow(gff))
  
  annot = data.frame(attr = strsplit(gff$attribute[1],";"),Key =gff$Key[1])
  colnames(annot) = c("Attr","Key")
  for( i in 2:nrow(gff))
  {
    tmp = data.frame(attr = strsplit(gff$attribute[i],";"),Key =gff$Key[i])
    colnames(tmp) = c("Attr","Key")
    annot = bind_rows(annot, tmp)
  }
  
  
  
  return(full_join(gff %>% select(-attribute), annot %>% separate(Attr,c("AttrType","AttrValue"), sep = "=") %>% spread(AttrType,AttrValue)) %>% select(-Key))
}