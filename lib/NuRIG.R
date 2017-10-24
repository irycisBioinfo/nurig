library(dplyr)
library(tidyr)
library(XML)

#setwd("~/ST117/Paper/data/ST117/")
 #listFile = "list2"
 #annotFile = "annot.tab"
 #ref = "CP018065.1EnterococcusfaeciumstrainE1,completegenome.fasta"
 #ref = "GCA_001562905.1_ASM156290v1_genomic.fna"
 outFile = "out.xml"

lista = read.table(listFile, sep = "\t")
#colnames(lista) = c("file")

if(exists("annotFile"))
{
  annot = read.delim(annotFile,sep = "\t", quote = "",stringsAsFactors = FALSE, header =FALSE)
}
if(!exists("palette_color"))
{
	palette_color = 'rainbow'
}

if (ncol(lista) == 2)
{
  colnames(lista) = c("file", "names")
} else{
  colnames(lista) = c("file")
  lista$names <- lista$file
}

print("Comparissons");

j = 1
for (i in lista$file)
{
  system(paste("nucmer ", ref, i, "> /dev/null", collapse = ""))
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

theTable = finalTable






finalTable = theTable

finalTable = finalTable %>% group_by(rID) %>% summarise(ln = first(lenref)) %>% mutate(prueba = cumsum(ln)) %>% mutate(position = prueba -ln) %>% select(rID, position) %>% full_join(finalTable, .)

legend = finalTable %>% group_by(file,names) %>% summarise(score = sum(peridentity * lenAlig /100)) %>% arrange(desc(score))
legend$order = seq(1:length(legend$file))

finalTable = finalTable %>% full_join(legend)

size = finalTable %>% group_by(rID) %>% summarise(tl = first(lenref)) %>% summarise(S = sum(tl))




xml = xmlTree()
xml$addTag(
  "cgview",
  attrs = c(
    backboneRadius = "160",
    sequenceLength = as.integer(size[1]),
    title = ref,
    labelPlacementQuality = "best"
  ),
  close = FALSE
)



j = 1

if (palette_color == 'rainbow')
{
  colores = rainbow(max(finalTable$slot))
} else if (palette_color == 'topo') {
  colores = topo.colors(max(finalTable$slot))
} else if (palette_color == 'terrain') {
  colores = terrain.colors(max(finalTable$slot))
} else if (palette_color == 'distinct'){
  colores = randomcoloR::distinctColorpalette_color(max(finalTable$slot))
} else {
  colores = colores = rainbow(max(finalTable$slot))
}

for (tab in 1:max(legend$order))
{
  tmp = finalTable %>% filter(order == tab)
  
  xml$addTag("featureSlot",
             attrs = c(strand = "direct"),
             close = FALSE)
  clr = paste("rgb(", paste(col2rgb(colores[j]), collapse = ","), ")", collapse = "")
  j = j + 1
  
  xml$addTag("feature",
             attrs = c(color = clr, decoration = "arc"),
             close = FALSE)
  for (r in 1:nrow(tmp))
  {
    #xml$addTag("feature", close = FALSE)
    xml$addTag(
      "featureRange",
      attrs = c(
        start = tmp$startreference[r] + tmp$position[r],
        stop = tmp$endreference[r] + tmp$position[r],
        opacity = tmp$peridentity[r] / 100
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


for (i in legend$names)
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

contigs = finalTable %>% group_by(rID) %>% summarise(position = first(position) +1)

lblBoolean = !(exists("annot"))

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
        showLabel = lblBoolean
      ),
      close = TRUE
    )
  }
  xml$addTag(
    "featureRange",
    attrs = c(
      start = contigs$position[nrow(contigs)],
      stop = as.integer(size[1]),
      mouseover = contigs$rID[nrow(contigs)],
      label = contigs$rID[nrow(contigs)],
      color = cl[i + 1],
      showLabel = lblBoolean
    ),
    close = TRUE
  )
  xml$closeTag()
  xml$closeTag()
}

### Draw Annotation

if(exists("annot"))
{
  if(ncol(annot) <6){
    annot$color = "black"
  }
  colnames(annot) = c("rID","start","end","strand","label","color")
  finalAnnot = finalTable %>% group_by(rID) %>% summarise(position = first(position)) %>% left_join(., annot, by ="rID") %>%  mutate(start = start +
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
	  
		  stop = tmp$end[i]
		  start = tmp$start[i]
	    if(stop > max(finalTable$lenref))
	    {
	      stop = max(finalTable$lenref)
	    }
		  if(start > max(finalTable$lenref))
		  {
		    start = max(finalTable$lenref)
		  }
    xml$addTag(
      "featureRange",
      attrs = c(
        start = start,
        stop = stop,
        color = tmp$color[i],
        label = tmp$label[i],
        showLabel = "TRUE"
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
    
    stop = tmp$end[i]
    start = tmp$start[i]
    if(stop > max(finalTable$lenref))
    {
      stop = max(finalTable$lenref)
    }
    if(start > max(finalTable$lenref))
    {
      start = max(finalTable$lenref)
    }
    xml$addTag(
      "featureRange",
      attrs = c(
        start = start,
        stop = stop,
        color = tmp$color[i],
        label = tmp$label[i],
        showLabel = "TRUE"
      ),
      close = TRUE
    )
  }
  xml$closeTag()
  xml$closeTag()
}
saveXML(xml$value(), file = "out.xml")
#write.table(finalTable,file = "FinalTable.tsv",sep = '\t', row.names = FALSE)

