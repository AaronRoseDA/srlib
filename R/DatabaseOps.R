
create_db_connection <- function(user_name, password) {
  con <- tryCatch({
    dbConnect(odbc::odbc(),
              Driver = "ODBC Driver 17 for SQL Server",
              Server = "symregserver.database.windows.net",
              Database = "SymReg",
              UID = user_name,
              PWD = password,
              Encrypt = "yes",
              TrustServerCertificate = "no")
  }, error = function(e) {
    cat("Error: Unable to connect to the database. Please check your credentials.\n")
    return(NULL)
  })

  if (!is.null(con)) {
    cat("Database connection established successfully.\n")
  }

  return(con)
}



library(odbc)

con <- dbConnect(odbc::odbc(),
                 Driver   = "ODBC Driver 18 for SQL Server",
                 Server   = "symregserver.database.windows.net",
                 Database = "SymReg",
                 UID      = "AaronRose",   # Not your Gmail
                 PWD      = "-2025-g!",
                 Encrypt  = "yes",
                 TrustServerCertificate = "no",
                 Timeout  = 30)

dbDisconnect(conn = con)

append_to_table <- function(connection, table_name, data, verbose = FALSE, DATE_INSERTED = TRUE, mode = "append") {

  # Validate mode
  if (!mode %in% c("append", "replace")) stop("Invalid mode. Use 'append' or 'replace'.")

  # Check if the table exists
  tables <- dbGetQuery(connection, "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA = 'dbo'")
  if (!(table_name %in% tables$TABLE_NAME)) {
    cat("Table not found. Available tables:\n")
    print(tables)
    table_name <- readline("Enter the correct table name (exact match required) or press Enter to cancel: ")
    if (!(table_name %in% tables$TABLE_NAME)) {
      cat("Operation canceled.\n")
      return(NULL)
    }
  }

  # Fetch the table structure
  db_columns <- dbListFields(connection, table_name)

  # Ensure DATE_INSERTED is in both input and database schema
  if (DATE_INSERTED) {
    if ("DATE_INSERTED" %in% db_columns) {
      data$DATE_INSERTED <- Sys.time()  # Update or add DATE_INSERTED
    } else {
      data$DATE_INSERTED <- NULL  # Remove it if not in DB schema
    }
  }

  # Ensure column names match
  data <- data[, intersect(names(data), db_columns), drop = FALSE]

  # Confirm before replacing data
  if (mode == "replace") {
    confirm <- readline(sprintf("Are you sure you want to replace all data in '%s'? Type 'yes' to confirm: ", table_name))
    if (tolower(confirm) != "yes") {
      cat("Operation canceled.\n")
      return(NULL)
    }
    dbExecute(connection, paste0("DELETE FROM " , table_name))
  }

  # Append the data
  dbWriteTable(connection, table_name, data, append = TRUE, row.names = FALSE)

  if (verbose) cat("Data successfully ", mode, "ed into table: ", table_name, "\n",sep = "")
}


create_table <- function(connection, table_name, table, DATE_INSERTED = TRUE) {
  # connection <-  con
  # table_name <- "input_NORMAL"
  # table<- NORMAL_data
  # DATE_INSERTED = FALSE
  # Validate input
  if (!is.data.frame(table)) stop("The 'table' parameter must be a data frame.")

  # Check if the table exists
  tables <- dbGetQuery(connection, "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA = 'dbo'")
  table_exists <- table_name %in% tables$TABLE_NAME

  if (table_exists) {
    # Table exists, prompt for append/replace
    confirm <- readline(sprintf("Table '%s' already exists. Append or replace? (append/replace/cancel): ", table_name))
    if (!confirm %in% c("append", "replace")) {
      cat("Operation canceled.\n")
      return(NULL)
    }
    append_to_table(connection, table_name, table, verbose = TRUE, DATE_INSERTED = DATE_INSERTED, mode = confirm)
  } else {
    # Create a new table
    # if (DATE_INSERTED && !"DATE_INSERTED" %in% names(table)) table$DATE_INSERTED <- Sys.time()
    if (nrow(table) > 0) {
      if (DATE_INSERTED && !"DATE_INSERTED" %in% names(table)) {
        table$DATE_INSERTED <- Sys.time()
      }
    } else {
      table$DATE_INSERTED <- character()
    }

    # Generate CREATE TABLE SQL
    create_stmt <- paste0("CREATE TABLE ", table_name, " (",
                          paste(sprintf("[%s] %s", names(table), ifelse(sapply(table, is.numeric), "FLOAT", "VARCHAR(255)")), collapse = ", "),
                          ");")

    # Execute table creation
    dbExecute(connection, create_stmt)

    # Insert data if available
    if (nrow(table) > 0) append_to_table(connection, table_name, table, verbose = TRUE, DATE_INSERTED = DATE_INSERTED, mode = "append")

    cat("Table '", table_name, "' created successfully.\n")
  }
}


delete_db_items <- function(connection, items = NULL, mode) {

  # Validate mode
  if (!mode %in% c("table", "view")) {
    stop("Invalid mode. Use 'table' or 'view'.")
  }

  # If no items are provided, get all objects of the specified mode
  if (is.null(items) || length(items) == 0) {
    query <- if (mode == "table") {
      "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA = 'dbo'"
    } else {
      "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.VIEWS WHERE TABLE_SCHEMA = 'dbo'"
    }
    items <- dbGetQuery(connection, query)$TABLE_NAME
  }

  # Iterate through each item
  for (item in items) {

    # Fetch first 10 rows before deletion
    preview_query <- sprintf("SELECT TOP 10 * FROM dbo.[%s]", item)
    preview_data <- tryCatch(dbGetQuery(connection, preview_query), error = function(e) NULL)

    if (!is.null(preview_data)) {
      print(preview_data)  # Show data preview
      cat(sprintf("\nAre you sure you want to delete '%s' (mode: %s)? Type 'yes' to confirm: ", item, mode))
      user_input <- readline()

      if (tolower(user_input) != "yes") {
        cat(sprintf("Skipping deletion of '%s'.\n", item))
        next  # Skip to the next item
      }
    } else {
      cat(sprintf("'%s' does not exist or is empty.\n", item))
      next
    }

    # If mode is table, drop foreign key constraints first
    if (mode == "table") {
      constraint_query <- sprintf("
        SELECT CONSTRAINT_NAME
        FROM INFORMATION_SCHEMA.TABLE_CONSTRAINTS
        WHERE CONSTRAINT_TYPE = 'FOREIGN KEY'
        AND TABLE_SCHEMA = 'dbo'
        AND TABLE_NAME = '%s'", item)

      constraints <- dbGetQuery(connection, constraint_query)

      if (nrow(constraints) > 0) {
        for (constraint in constraints$CONSTRAINT_NAME) {
          drop_constraint_query <- sprintf("ALTER TABLE dbo.[%s] DROP CONSTRAINT [%s];", item, constraint)
          try(dbExecute(connection, drop_constraint_query), silent = TRUE)
        }
      }

      # Drop the table
      delete_query <- sprintf("DROP TABLE dbo.[%s];", item)
    } else {  # Mode is 'view'
      delete_query <- sprintf("DROP VIEW dbo.[%s];", item)
    }

    # Execute deletion
    tryCatch({
      dbExecute(connection, delete_query)
      cat(sprintf("'%s' deleted successfully.\n", item))
    }, error = function(e) {
      cat(sprintf("Error deleting '%s': %s\n", item, e$message))
    })
  }

  cat("Deletion process completed.\n")
}


get_db_table <- function(connection, table = NULL, n_rows = NULL, summary = FALSE, verbose = FALSE) {

  # If summary = TRUE, ignore table & n_rows and return all tables with first 10 rows
  if (summary) {
    tables <- dbGetQuery(connection,
                         "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA = 'dbo'")
    table_list <- list()


    if (is.null(table))
      for (tbl in tables$TABLE_NAME) {
        query <- sprintf("SELECT TOP 4 * FROM dbo.[%s]", tbl)
        table_list[[tbl]] <- dbGetQuery(connection, query)
      } else if (table %in% tables$TABLE_NAME){
        query <- sprintf("SELECT TOP 4 * FROM dbo.[%s]", table)
        table_list[[table]] <- dbGetQuery(connection, query)
      }

    if (verbose){
      print(table_list)}
    return(table_list)
  }

  # Validate table parameter
  if (is.null(table)) {
    stop("Please specify a valid table name or set summary = TRUE.")
  }

  # Check if table exists
  tables <- dbGetQuery(connection,
                       "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA = 'dbo'")

  if (!(table %in% tables$TABLE_NAME)) {
    cat("Invalid table name. Available tables:\n")
    print(tables)
    return(NULL)
  }

  # Construct query based on n_rows parameter
  query <- if (is.null(n_rows)) {
    sprintf("SELECT * FROM dbo.[%s]", table)
  } else {
    sprintf("SELECT TOP %d * FROM dbo.[%s]", n_rows, table)
  }

  # Fetch data
  data <- dbGetQuery(connection, query)
  if (verbose){
    print(head(data))
  }

  # If n_rows is specified, print table structure and head
  if (!is.null(n_rows)) {
    cat("\nTable Structure:\n")
    str(data)
    cat("\nTable Head:\n")
    print(head(data))
  }

  return(data)
}





delete_db_items(con,items = "input_NORMAL",mode = "table")
delete_db_items(con,items = "input_UNIFORM",mode = "table")
delete_db_items(con,items = "input_TRIANGULAR",mode = "table")
create_table(con, "input_TRIANGULAR", tri, DATE_INSERTED = FALSE)
nor[sapply(nor, is.infinite)] <- NA
create_table(con, "input_NORMAL", nor, DATE_INSERTED = FALSE)
create_table(con, "input_UNIFORM", uni, DATE_INSERTED = FALSE)






