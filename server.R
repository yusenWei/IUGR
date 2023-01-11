#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(Seurat)
library(ggpubr)
options(shiny.maxRequestSize=200000*1024^2)
scRNA<- readRDS("data/scRNA.rds")
scRNA1<- readRDS("data/human_pig.rds")
DefaultAssay(scRNA) <- "RNA"
scRNA@assays$RNA@data <- scRNA@assays$RNA@counts
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
DefaultAssay(scRNA1) <- "RNA"
scRNA1@assays$RNA@data <- scRNA1@assays$RNA@counts
scRNA1 <- NormalizeData(scRNA1, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA1 <- ScaleData(scRNA1, features = rownames(scRNA1))

# Define server logic required to draw a histogram
color_cluster=c("#6AB3F4","#D459A2","#3A8F93","#FACC13",
"#5A7FF8","#3A3AC8",
"#8747E0","#12C2C4","#FF7B78",
"#B27EEC","#CE1222",
"#45CFCF","#EBA8A6","#FF6347","#F14764",
"#159CA1","#A21A6D",
"#EBA8A6","#e9148f",
"#AE6800","#32C35E","#5AC722")
names(color_cluster)=c("NK/NKT cells", "T cells","Proliferation","Cycling NKT cells",
"Neutrophil-myeloid progenitors","Cycling T cells",
"Monocytes", "cDC2","cDC1",
"B cells","Plasma B cells",
"Mesenchymal cells","HSC","Endothelial cells","Kupffer cells",
"Neutrophils","Myelocytes",
"Hepatocytes","Cholangiocytes",
"pDC","Basophils","Erythroid cells")
color_cluster1=c("#6AB3F4","#D459A2",
"#5A7FF8",
"#8747E0","#12C2C4","5AC722",
"#B27EEC",
"#45CFCF","#EBA8A6","#FF6347","#F14764",
"#159CA1","#A21A6D",
"#EBA8A6","#e9148f",
"#32C35E","#5AC722")
names(color_cluster1)=c("NK/NKT cells", "T cells",
"Neutrophil-myeloid progenitors",
"Monocytes", "cDC","Late erythroid",
"B/Plasma cells",
"Mesenchymal cells","HSC","Endothelial cells","Macrophage cells",
"Neutrophils","Myelocytes",
"Hepatocytes","Cholangiocytes",
"Basophils","Early/Mid erythroid")
color_group=c("#98BCDA","#EBA8A6")
names(color_group)=c("NBW","IUGR")
resubset<-function(scRNA,gender,cluster){
  Idents(scRNA) <- "gender"
  scRNA <- subset(scRNA,idents = gender)
  Idents(scRNA) <- "celltype"
  scRNA <- subset(scRNA,idents = cluster)
}
plot<-function(scRNA,gene,plot){
  if(plot == "Vlnplot"){
    VlnPlot(scRNA,features =gene,split.by = "group",pt.size = 0.01)+theme_bw()+
      theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))}
  else if(plot == "FeaturePlot"){
    FeaturePlot(scRNA,features =gene,split.by = "group",pt.size = 2)}
  else if(plot == "Dotplot"){
    DotPlot(scRNA,features =gene,split.by = "group")+theme_bw()+coord_flip()}
}


shinyServer(function(input, output){
    # Create the UMAP visualization and render it in the plot
    output$umap <- renderPlot({
      DimPlot(scRNA, reduction = "umap",label = T,label.box = T,cols = color_cluster,repel = T)+ NoLegend()
    })
    doplot <- eventReactive(input$action,{
      plot(resubset(scRNA,input$gender,input$cluster),input$gene,input$method)})
    output$plot1 <- renderPlot({doplot()})
    output$cluster<-renderUI({
      clusterlist<-unique(scRNA@meta.data$celltype)
      selectInput("cluster",'Please select the cell type you are interested in:',list=as.vector(clusterlist))
    })
    output$gender<-renderUI({
      genderlist<-unique(scRNA@meta.data$gender)
      selectInput("gender",'Please select the gender you are interested in:',list=as.vector(genderlist))
    })
    output$method<-renderUI({
      selectInput("method",'Please select the method to display',list=c("Vlnplot","FeaturePlot","Dotplot"))
    })
    output$umap1<- renderPlot({
      DimPlot(scRNA1, reduction = "umap",label = T,label.box = T,cols = color_cluster1,repel = T)+ NoLegend()
    })
    doplot1 <- eventReactive(input$action1,{
      plot(subset(scRNA1,idents=input$cluster1),input$gene1,input$method1)})
    output$plot2 <- renderPlot({doplot1()})
    output$dataset <- DT::renderDataTable({
      markers<- read.csv("data/markers.csv",row.names = 1)
      DT::datatable(markers,extensions = 'Buttons',
                    options = list(dom='Bfrtip',
                                   buttons=c('copy', 'csv', 'excel', 'print', 'pdf')),
                    caption="Table 1. irisdata",filter = "top")})
    output$dataset1 <- DT::renderDataTable({
      markers1<- read.csv("data/markers1.csv",row.names = 1)
      DT::datatable(markers1,extensions = 'Buttons',
                    options = list(dom='Bfrtip',
                                   buttons=c('copy', 'csv', 'excel', 'print', 'pdf')),
                    caption="Table 2. irisdata",filter = "top")})
    output$plotDown <- downloadHandler(
      filename = function(){
        paste0(input$gene, '.',input$extPlot)
      },
      content = function(file){
        if(input$extPlot == 'pdf'){
          pdf(file)
        }else if(input$extPlot == 'png'){
          png(file)
        }else{
          jpeg(file)
        }
        print(plot(resubset(scRNA,input$gender,input$cluster),input$gene,input$method))
        dev.off()
      }
    )
    output$plotDown1 <- downloadHandler(
      filename = function(){
        paste0(input$gene1, '.',input$extPlot1)
      },
      content = function(file){
        if(input$extPlot == 'pdf'){
          pdf(file)
        }else if(input$extPlot == 'png'){
          png(file)
        }else{
          jpeg(file)
        }
        print(plot(subset(scRNA1,idents=input$cluster1),input$gene1,input$method1))
        dev.off()
      }
    )
  })
