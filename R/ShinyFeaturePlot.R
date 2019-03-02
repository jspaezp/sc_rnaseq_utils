library(shiny)
require(Seurat)
library(DT)



#' Gets a data frame from a seurat object
#' 
#' Provided a seurat object, returs a data frame of the count values, being 
#' the columns each 'gene' and the rows each UMI/cell.
#'
#' It returns only the genes annotated as variable and the identity column.
#'
#' @param seurat A seurat object
#'
#' @return
#' @export
#'
#' @examples
#' > as.data.frame.Seurat(pbmc_small)[1:3,1:3]
#'                    LTB EAF2 CD19
#' ATGCCAGAACGACT 6.062788    0    0
#' CATGGCCTGTGCAT 6.714813    0    0
#' GAACCTGATGAACC 7.143118    0    0
#' @importFrom Seurat FetchData GetIdent
as.data.frame.Seurat <- function(seurat, vars.all = seurat@var.genes) {
    tmp <- Seurat::FetchData(seurat, vars.all = vars.all) 
    tmp <-  as.data.frame(data.matrix(tmp), stringsAsFactors = FALSE) 
    
    tmp$ident <- Seurat::GetIdent(seurat, uniq = FALSE, cells.use = rownames(tmp))
    return(tmp)
}


#' Simulates flow cytometry plotting
#' 
#' Note that it imputes a little bit of noise to values of 0 just for the
#' Sake of visualization
#'
#' @param object data
#' @param markernames 
#'
#' @return a ggplot grid with the plots
#' @export
#'
#' @examples
plot.flowstyle <- function(object, markernames) {
    UseMethod("plot.flowstyle", object)
} 


#' @importFrom GGally ggpairs wrap
plot_flowstyle.data.frame <- function(object, markernames, classif_col = "ident") {

    object <- object[,c(markernames, classif_col)] 
    tmp_ident <- object$ident
    
    suppressWarnings({
        object[object == 0] <- abs(runif(sum(object == 0),min = 0, max = 0.4))
        object$ident <- tmp_ident
    })
    
    g <- GGally::ggpairs(
        as.data.frame(object),
        columns = 1:(ncol(object) - 1),
        ggplot2::aes_string(colour = "ident"),
        progress = FALSE,
        lower = list(
            continuous = GGally::wrap(
                "dot_no_facet",
                alpha = 0.2)),
        diag = list(
            continuous = GGally::wrap(
                'densityDiag',
                alpha = 0.3)),
        upper = list(
            continuous =  GGally::wrap(
                "density",
                alpha = 0.4))) + 
        ggplot2::theme_bw()
    print(g)
}



plot.flowstyle.seurat <- function(object, markernames, classif_col = "ident") {
    tmp <- as.data.frame.Seurat(object, markernames)
    plot_flowstyle.data.frame(tmp, markernames = markernames,
                              classif_col = classif_col)
}


ranger.seurat <- function(object, pval_cutoff = 0.05, imp_method = c("janitza", "altmann"), ...) {
    
    base_ranger <- function(importance, ...) {
        ranger_fit <- ranger::ranger(
            ident ~ .,
            data = tmp,
            num.trees = 1000, 
            mtry = floor(ncol(tmp)/5),
            importance = "permutation",
            ...)
    }
    
    tmp <- as.data.frame.Seurat(object, object@var.genes)
    
    if (imp_method == "altmann") {
        ranger_fit <- base_ranger(importance = imp_method, ...)
        importances_ranger <- ranger::importance_pvalues(
            ranger_fit, "altman", 
            formula = ident ~ ., data = tmp)
        
    } else if (imp_method == "janitza") {
        ranger_fit <- base_ranger(importance = imp_method, ...)
        importances_ranger <- ranger::importance_pvalues(ranger_fit)
        
    } else {
        stop("'imp_method' does not match either 'altmann' or 'janitza'")
    }
    
    signif_importances_ranger <- importances_ranger[
        importances_ranger[,'pvalue'] < pval_cutoff,]
    
    signif_importances_ranger <- signif_importances_ranger[
        order(signif_importances_ranger[,"importance"], decreasing = TRUE), ]
    
    signif_importances_ranger <- as.data.frame(signif_importances_ranger)
    signif_importances_ranger[["gene"]] <- rownames(signif_importances_ranger)
    
    return(list(ranger_fit, importances_ranger, signif_importances_ranger))
}


#' Title
#'
#' @param seurat_object 
#' @param starting_genes 
#'
#' @return
#' @export
#'
#' @examples
featureplot_gadget <- function(seurat_object = pbmc_small,
                               starting_genes = NULL) {
    
    plotting_pannel <- sidebarPanel(
        h1("Plotting single cell data !"),
        uiOutput("allgenes"),
        actionButton("add_gene", "Add to plottable genes"),
        uiOutput("genes"),
        actionButton("plot_command", "PLOT!"),
        checkboxInput(
            "hidelegends",
            "Hide Legends",
            value = FALSE),
        sliderInput(
            "pointsize",
            "Point Size:",
            min = 0.5,
            max = 10,
            value = 4
        )
    )
    
    
    ui <- shiny::fluidPage(
        # Application title
        shiny::titlePanel("Feature Plot for a Seurat Object"),

        shiny::sidebarLayout(
            plotting_pannel,

            # Show a plot of the generated distribution
            shiny::mainPanel(
                shiny::tabsetPanel(
                    shiny::tabPanel(
                        "Find Markers", 
                        shiny::div(DT::dataTableOutput("rangerImpTable")),
                        shiny::plotOutput("rangerImpPlot")),
                    shiny::tabPanel("TSNE plot", shiny::plotOutput("tsnePlot", height = '600px')),
                    shiny::tabPanel("Violin Plot", shiny::plotOutput("violinPlot", height = '600px')),
                    shiny::tabPanel("Pairs Plot", shiny::plotOutput("pairsPlot", height = '800px')))
            )
            
        )
    )

    server <- function(input, output) {

        colorscale <- viridis::viridis(2, direction = -1)
        total_genes <- rownames(seurat_object@raw.data)
        
        importance_df <- ranger.seurat(seurat_object)[[3]]
       
        if (is.null(starting_genes)) {
            plottable_genes <- rownames(head(importance_df))
        } else {
            plottable_genes <- starting_genes
        }

        get_plottable_genes <- shiny::eventReactive(input$add_gene, {
            plottable_genes <<- sort(unique(c(
                        plottable_genes, input$allgenes)))
            return(plottable_genes)
        })
        
        output$rangerImpTable <- DT::renderDataTable({
            importance_df
        })
        
        output$rangerImpPlot <- shiny::renderPlot({
            importance_df[["gene"]] <-  forcats::fct_reorder(
                importance_df[["gene"]], 
                importance_df[["importance"]])
            
            g <- ggplot2::ggplot(
                importance_df, 
                ggplot2::aes_string(x = "importance", xend = 0,
                           y = "gene", yend = "gene")) + 
                ggplot2::geom_point() + 
                ggplot2::geom_segment() + 
                ggplot2::ggtitle("Relative variable importance")
            g
        })

        output$genes <- shiny::renderUI({
            shiny::checkboxGroupInput(
                inputId = "genenames",
                label = "Genes to plot",
                choices = get_plottable_genes(),
                selected = get_plottable_genes())
        })

        output$allgenes <- shiny::renderUI({
            shiny::selectInput(
                'allgenes',
                "allgenes",
                choices = sort(
                    total_genes[!total_genes %in% get_plottable_genes()]))
        })

        captured_input <- shiny::eventReactive(input$plot_command, {
            retlist <- list(
               genenames = input$genenames,
               numcolors = input$numcolors,
               pointsize = input$pointsize,
               hidelegends = input$hidelegends)
            return(retlist)
        })

       output$tsnePlot <- shiny::renderCachedPlot({
           capt_input <- captured_input()

           Seurat::FeaturePlot(
                object = seurat_object,
                features.plot = capt_input$genenames,
                cols.use = colorscale,
                reduction.use = "tsne",
                dim.1 = 1,
                no.legend = capt_input$hidelegends,
                do.return = FALSE,
                pt.size = capt_input$pointsize)
       }, cacheKeyExpr = {
           captured_input()
       })
       
       output$violinPlot <- shiny::renderCachedPlot({
           capt_input <- captured_input()
           
           Seurat::VlnPlot(
               object = seurat_object,
               features.plot = capt_input$genenames,
               do.return = FALSE)
       }, cacheKeyExpr = {
           captured_input()
       })
       
       output$pairsPlot <- shiny::renderCachedPlot({
           capt_input <- captured_input()
           
           plot.flowstyle(
               object = seurat_object,
               markernames = make.names(capt_input$genenames))
       }, cacheKeyExpr = {
           captured_input()
       })
       
    }

    # Run the application
    shiny::shinyApp(ui = ui, server = server)

}

featureplot_gadget()

