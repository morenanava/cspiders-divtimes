###############################################################################
# The below is a custom taxon range through time plot, modified after the     #
# function within the Palaeoverse. This allows a range plot ordered by groups #
# above species. Added  title = "Temporal range of taxa", plot_names = FALSE, #
# group_by = NULL, x_pad = 0                                                  #
###############################################################################
tax_range_time_RJG  <- function( occdf, name = "genus", min_ma = "min_ma",
                                 max_ma = "max_ma", by = "FAD", plot = FALSE,
                                 title = "Temporal range of taxa",
                                 plot_names = FALSE, group_by = NULL,
                                 x_pad = 0, x_range = NULL, y_range = NULL,
                                 margin_text = TRUE, xlabline = 4,
                                 ital_font = TRUE )
{
  if( is.data.frame( occdf ) == FALSE ){
    stop( "`occdf` should be a dataframe" )
  }
  if( is.logical( plot ) == FALSE ){
    stop( "`plot` should be logical (TRUE/FALSE)" )
  }
  if( ! is.numeric( occdf[, max_ma] ) || ! is.numeric( occdf[, min_ma] ) ){
    stop( "`max_ma` and `min_ma` must be of class numeric." )
  }
  if( any( c( name, min_ma, max_ma, group_by ) %in% colnames( occdf ) == FALSE ) ){
    stop( "Either `name`, `min_ma`, `max_ma`, or `group_by`\n",
          "is not a named column in\n`occdf`" )
  }
  if( any( is.na( occdf[, name] ) ) ){
    stop( "The `name` column contains NA values" )
  }
  if( any( is.na( occdf[, min_ma] ) ) || any( is.na( occdf[, max_ma] ) ) ){
    stop( "`min_ma` and/or `max_ma` columns contain NA values" )
  }
  if( ! by %in% c( "name", "FAD", "LAD" ) ){
    stop( "`by` must be either \"FAD\", \"LAD\", or \"name\"" )
  }
  if( any( is.na(occdf[, group_by])) ){
    stop( "`group_by` contains NA values" )
  }
  unique_taxa <- unique( occdf[, name] )
  unique_taxa <- unique_taxa[order( unique_taxa )]
  temp_df <- data.frame( taxon = unique_taxa,
                         taxon_id = seq( 1, length( unique_taxa ), 1 ),
                         max_ma = rep( NA, length( unique_taxa ) ),
                         min_ma = rep( NA, length( unique_taxa ) ),
                         range_myr = rep( NA, length( unique_taxa ) ),
                         n_occ = rep( NA, length( unique_taxa ) ),
                         group_by = rep( NA, length( unique_taxa ) ) )
  for( i in seq_along( unique_taxa ) ){
    vec <- which( occdf[, name] == unique_taxa[i] )
    temp_df$max_ma[i] <- max( occdf[vec, max_ma] )
    temp_df$min_ma[i] <- min( occdf[vec, min_ma] )
    temp_df$range_myr[i] <- temp_df$max_ma[i] - temp_df$min_ma[i]
    temp_df$n_occ[i] <- length( vec )
    #RJG - If we need to group by something, we need to add this to DF
    if( ! is.null( group_by ) ) temp_df$group_by[i] <- occdf[vec[1], group_by]
  }
  row.names( temp_df ) <- NULL
  str2round <-  temp_df[,c( "max_ma","min_ma","range_myr" )]
  temp_df[, c( "max_ma", "min_ma", "range_myr" )] <- round( x = str2round,
                                                            digits = 3 )
  if( by == "FAD" ){
    temp_df <- temp_df[order( temp_df$max_ma ), ]
    temp_df$taxon_id <- seq_len( nrow( temp_df ) )
  }
  if( by == "LAD" ) {
    temp_df <- temp_df[order( temp_df$min_ma ), ]
    temp_df$taxon_id <- seq_len( nrow( temp_df ) )
  }
  #RJG - now if we need to group by: since this is already ordered by FAD or
  #LAD that is carried through as the order within each group_by category
  if( ! is.null( group_by ) ){
    temp_df <- temp_df[order( temp_df$group_by ), ]
    temp_df$taxon_id <- seq_len( nrow( temp_df ) )
  }

  if( plot == TRUE ){
    if( is.null( x_range ) ) x_range <- c( max( temp_df$max_ma ),
                                           min( temp_df$min_ma ) )
    x_range[2] <- x_range[2] - x_pad
    if( is.null( y_range ) ) y_range <- c (0, nrow( temp_df ) )

    #RJG - If plotting names, remove Y axis and labels, then add text using
    #text. Position this to right of min age. For lots of taxa needs to be
    #output to a device defining a sensible aspect ratio
    if( plot_names == TRUE ){

      #RJG - If we are adding names in the margin, we can expand outer margin
      #area to allow these, but we need to do this before we call plot
      if( margin_text == TRUE && ! is.null( group_by ) ){
        maxLength <- max( nchar( temp_df$group_by ) )
        #This is about the right amount unless modifying with cex
        par( oma = c( 0,0,0,maxLength/3 ) )
      }

      plot( x = NA, y = NA, xlim = x_range, ylim = y_range,axes = TRUE,
            xaxt = "n", xlab = NA, yaxt = "n", ylab = NA, main = title )

      #RJG - add grey boxes
      periods <- GTS2020[which( GTS2020$rank=="period" ),]
      for( i in 1:length( periods ) ){
        if( i%%2 == 0 ){
          graphics::rect( periods[i,]$max_ma, graphics::par()$usr[3],
                          periods[i,]$min_ma, graphics::par()$usr[4],
                          col = c( "grey90", "white" ), border = NA )
        }
      }
      #Add a border over the grey boxes to make the axis look OK
      graphics::rect( graphics::par()$usr[1], graphics::par()$usr[3],
                      graphics::par()$usr[2], graphics::par()$usr[4],
                      col = NA, border = "black" )
      
      #SAC: plotting text in italics
      if( ital_font == TRUE ){
        # Source: https://stackoverflow.com/questions/29943251/displaying-values-from-a-character-vector-as-italic-labels-in-boxplot-in-r
        make.italic <- function( x ){
          as.expression( lapply( x, function( y ) bquote( italic( .(y) ) ) ) )
        }
        itvec <- make.italic( temp_df$taxon )
        text( x = temp_df$min_ma, y = temp_df$taxon_id,
              itvec,
              cex = 0.8, pos = 4 )
      }else{
        text( x = temp_df$min_ma, y = temp_df$taxon_id, temp_df$taxon,
              cex = 0.8, pos = 4 )
      }
      
      # Add labels for group_by
      if( ! is.null( group_by ) ){
        unique_groups <- unique( temp_df$group_by )
        for( i in seq_along( unique_groups ) ){
          vec <- which( temp_df$group_by == unique_groups[i] )
          maxLength <- max( nchar( temp_df$taxon[vec] ) ) + 2
          minMa <- min( temp_df$min_ma[vec] )
          max <- max( vec )
          min <- min( vec )
          if( margin_text == TRUE ) mtext( unique_groups[i], side = 4,
                                           at = (min + max)/2, las = 1,
                                           line = 0.1 )
          else text( x = minMa - (maxLength*3), y = (min + max)/2,
                     unique_groups[i], pos = 4 )
          if( min > 1 ) abline( h = min -.5, col="grey" )
        }
      }
    }
    else{
      plot( x = NA, y = NA, xlim = x_range, ylim = y_range, axes = TRUE,
            xaxt = "n", xlab = NA, ylab = "Taxon ID", main = title )
    }
    segments( x0 = temp_df$max_ma, x1 = temp_df$min_ma, y0 = temp_df$taxon_id,
              col = "black" )
    points( x = temp_df$max_ma, y = temp_df$taxon_id, pch = 20, col = "black" )
    points( x = temp_df$min_ma, y = temp_df$taxon_id, pch = 20, col = "black" )

    #RJG - customisations to plot for his publication - change dates to
    #intervals, do 1DP, and move closer to axis
    axis_geo( side = 1, intervals = "periods", exact = TRUE, round = 1,
              cex.axis = 0.9, padj = -1 )
    title( xlab = "Time (Ma)", line = xlabline )
    
  }
  return( temp_df )
}