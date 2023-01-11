#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
scRNA<- readRDS("data/scRNA.rds")
scRNA1<- readRDS("data/human_pig.rds")

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  # Add a plot of the UMAP visualization
    navbarPage('WangLab',inverse = T,collapsible = T,
             tabPanel("Guide",titlePanel("Comprehensive single-cell transcriptional profile of IUGR piglet liver with high homology to human")
                      ,h2("Vignettes"),h3("Intrauterine growth restriction (IUGR) is a problem in both human medicine and animal husbandry, and it has a direct impact on both the health and economy. This urgent event requires us to uncover the molecular mechanisms of the disease in order to create more potent interventional strategies. Herein, we constructed single-cell transcriptome profiles (n=41,969) of IUGRs and NBWs (female n=1, male n=3) one-week after birth. We discussed the difference in male and female IUGR, and established that pigs could sever as a useful model for researching on IUGR in humans.
                                          "),h3("On this site, you may explore the variation in expression of IUGRs and NBWs across various subpopulations on the 'IUGR Piglets' page. Additionally on the 'Human and pigs' page, where the list of genes could be retrieved and downloaded, you could examine the levels of gene expression within different subpopulations of pigs and humans.")
                      ,h2("Statement"),h3("The dataset combines the following two studies."),
                      h4("1.Popescu DM, Botting RA, Stephenson E, Green K, Webb S, Jardine L, Calderbank EF, Polanski K, Goh I, Efremova M, Acres M, Maunder D, Vegh P, Gitton Y, Park JE, Vento-Tormo R, Miao Z, Dixon D, Rowell R, McDonald D, Fletcher J, Poyner E, Reynolds G, Mather M, Moldovan C, Mamanova L, Greig F, Young MD, Meyer KB, Lisgo S, Bacardit J, Fuller A, Millar B, Innes B, Lindsay S, Stubbington MJT, Kowalczyk MS, Li B, Ashenberg O, Tabaka M, Dionne D, Tickle TL, Slyper M, Rozenblatt-Rosen O, Filby A, Carey P, Villani AC, Roy A, Regev A, Roberts I, Behjati S, Laurenti E, Teichmann SA, Haniffa M. Decoding human fetal liver haematopoiesis. Nature. 2019 Oct;574(7778):365-371. doi: 10.1038/s41586-019-1652-y. Epub 2019 Oct 9. PMID: 31597962; PMCID: PMC6861135.
                          "),h4("2.Ramachandran P, Dobie R, Wilson-Kanamori JR, Dora EF, Henderson BEP, Luu NT, Portman JR, Matchett KP, Brice M, Marwick JA, Taylor RS, Efremova M, Vento-Tormo R, Carragher NO, Kendall TJ, Fallowfield JA, Harrison EM, Mole DJ, Wigmore SJ, Newsome PN, Weston CJ, Iredale JP, Tacke F, Pollard JW, Ponting CP, Marioni JC, Teichmann SA, Henderson NC. Resolving the fibrotic niche of human liver cirrhosis at single-cell level. Nature. 2019 Nov;575(7783):512-518. doi: 10.1038/s41586-019-1631-3. Epub 2019 Oct 9. PMID: 31597160; PMCID: PMC6876711.")
                      ), 
             tabPanel("IUGR Piglets",plotOutput("umap",width='60%',height='800px'),
                      sidebarLayout(
                        sidebarPanel(
                          conditionalPanel(condition = "input == ture",
                                           selectInput("gene", "Please select the genes you are interested in:",
                                                       c(as.list(rownames(scRNA@assays$RNA@counts)),"nCount_RNA","nFeature_RNA",
                                                         "percent_mito","percent_ribo","percent_hb"))),
                          conditionalPanel(condition = "input.cluster==true",
                                           selectInput("cluster", "Please select the cell type you are interested in:",
                                                       as.list(levels(scRNA)))),
                          conditionalPanel(condition = "input.gender==true",
                                           selectInput("gender", "Please select the gender you are interested in:",
                                                       choices = c("female","male"))),
                          conditionalPanel(condition = "input.method==true",
                                           selectInput("method", "Please select the method to display",
                                                       choices = c("Vlnplot","FeaturePlot","Dotplot"))),
                          actionButton("action", "Plot Expression",icon=icon('angle-double-right')),
                          radioButtons('extPlot', 'Plot output format',
                                       choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                          helpText('Please choose the format of table and plot that you need, the download 
               buttons are placed in respective tabs'),
                          tags$style('#plotDown {background-color: red; color: white}')
                        ),
                        mainPanel(
                          plotOutput("plot1"),downloadButton('plotDown',label="Download Plot"),
                          DT::dataTableOutput("dataset")))),
             
             tabPanel("Human and pigs",plotOutput("umap1",width='60%',height='800px'),
                      sidebarLayout(
                        sidebarPanel(
                          conditionalPanel(condition = "input == ture",
                                           selectInput("gene1", "Please select the genes you are interested in:",
                                                       c(as.list(rownames(scRNA1@assays$RNA@counts)),"nCount_RNA","nFeature_RNA",
                                                         "percent_mito","percent_ribo","percent_hb"))),
                          conditionalPanel(condition = "input.cluster==true",
                                           selectInput("cluster1", "Please select the cell type you are interested in:",
                                                       as.list(levels(scRNA1)))),
                          conditionalPanel(condition = "input.method==true",
                                           selectInput("method1", "Please select the method to display",
                                                       choices = c("Vlnplot","FeaturePlot","Dotplot"))),
                          actionButton("action1", "Plot Expression",icon=icon('angle-double-right')),
                          radioButtons('extPlot1', 'Plot output format',
                                       choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                          helpText('Please choose the format of table and plot that you need, the download 
               buttons are placed in respective tabs'),
                          tags$style('#plotDown1 {background-color: red; color: white}')
                        ),
                        mainPanel(
                          plotOutput("plot2"),downloadButton('plotDown1',label="Download Plot"),
                          DT::dataTableOutput("dataset1")))))
))
