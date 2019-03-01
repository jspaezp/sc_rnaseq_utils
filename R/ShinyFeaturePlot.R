library(shiny)
require(Seurat)
require(RColorBrewer)


plotting_pannel <- sidebarPanel(
    h1("PLOT !"),
    actionButton("plot_command", "PLOT!"),
    h2("Functional Settings"),
    uiOutput("allgenes"),
    actionButton("add_gene", "Add to plottable genes"),
    uiOutput("genes"),
    h2("Aesthetic Settings"),
    h3("Colors"),
    selectInput("palette_type",
                "Palette Type",
                choices = c(
                    'Divergent' = 'div',
                    'Sequential' = 'seq',
                    'Categorical' = 'qual'),
                selected = 'seq'),
    sliderInput(
        "numcolors",
        "Number of colors:",
        min = 2,
        max = 10,
        value = 5
    ),
    uiOutput("palettes"),
    h3("Etc"),
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


featureplot_gadget <- function(seurat_object = pbmc_small, starting_genes = NULL) {
    palettes <- RColorBrewer::brewer.pal.info[
        RColorBrewer::brewer.pal.info$colorblind,]

    ui <- fluidPage(
        # Application title
        titlePanel("Feature Plot for a Seurat Object"),

        sidebarLayout(
            plotting_pannel,

            # Show a plot of the generated distribution
            mainPanel(
                tabsetPanel(
                tabPanel("tab 1", plotOutput("tsnePlot", height = '600px')),
                tabPanel("tab 2", plotOutput("violinPlot", height = '600px')),
                tabPanel("tab 3", "contents"))
            )
            
        )
    )

    server <- function(input, output) {

        output$palettes <- renderUI({
            pal <- palettes[palettes$category == input$palette_type,]
            selectInput(
                "palette",
                "Color Palette:",
                choices = rownames(pal)
            )
        })

        total_genes <- rownames(seurat_object@raw.data)
        plottable_genes <- ifelse(
            is.null(starting_genes), 
            yes = c(sample(total_genes, 2)),
            no = starting_genes)

        get_plottable_genes <- eventReactive(input$add_gene, {
            plottable_genes <<- sort(unique(c(
                        plottable_genes, input$allgenes)))
            return(plottable_genes)
        })

        output$genes <- renderUI({
            checkboxGroupInput(
                inputId = "genenames",
                label = "Genes to plot",
                choices = get_plottable_genes(),
                selected = get_plottable_genes())
        })

        output$allgenes <- renderUI({
            selectInput(
                'allgenes',
                "allgenes",
                choices = sort(
                    total_genes[!total_genes %in% get_plottable_genes()]))
        })

        captured_input <- eventReactive(input$plot_command, {
            retlist <- list(
               genenames = input$genenames,
               numcolors = input$numcolors,
               palette = input$palette,
               pointsize = input$pointsize,
               hidelegends = input$hidelegends)
            return(retlist)
        })

       output$tsnePlot <- renderCachedPlot({
           capt_input <- captured_input()

            colorscale <- RColorBrewer::brewer.pal(
                                capt_input$numcolors,
                                capt_input$palette)

            FeaturePlot(
                object = pbmc_small,
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
       
       output$violinPlot <- renderCachedPlot({
           capt_input <- captured_input()
           
           Seurat::VlnPlot(
               object = pbmc_small,
               features.plot = capt_input$genenames,
               do.return = FALSE)
       }, cacheKeyExpr = {
           captured_input()
       })
    }

    # Run the application
    shinyApp(ui = ui, server = server)

}

featureplot_gadget()

