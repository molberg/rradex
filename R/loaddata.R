loaddata <- function(molfile) {
    lines <- readLines(molfile, n=-1)
    comments <- grep("!", lines)
    lines <- lines[-comments]

    molecule <- lines[1]
    amass <- as.numeric(lines[2])

    current <- 3
    nlev <- as.integer(lines[current])
    ## print(nlev)
    con <- textConnection(lines[seq(nlev)+current])
    levels <- read.delim(con, sep="", header=FALSE)
    close(con)
    levels <- levels[,-1]
    if (ncol(levels) > 3) {
        m <- ncol(levels)
        qnum <- do.call(paste, c(levels[,3:m], sep=""))
        levels <- cbind(levels[,1:2], qnum)
    }
    names(levels) <- c("energy", "weight", "qnum")
    levels$qnum <- as.factor(levels$qnum)

    current <- current+nlev+1
    ntrans <- as.integer(lines[current])
    ## print(ntrans)
    con <- textConnection(lines[seq(ntrans)+current])
    transitions <- read.delim(con, sep="", header=FALSE)
    close(con)
    transitions <- transitions[,-1]
    names(transitions) <- c("up", "low", "E.A","GHz","E.up")

    current <- current+ntrans+1
    npart <- as.integer(lines[current])
    ## print(npart)

    collisions <- list()
    for (ipart in seq(npart)) {
        current <- current+1
        kind <- as.integer(unlist(strsplit(lines[current]," "))[1])
        ## print(kind)
        partner <- c("H2", "o-H2", "p-H2", "e", "H", "He", "H+")[kind]
        ## print(partner)

        current <- current+1
        ncoll <- as.integer(lines[current])
        ## print(ncoll)

        current <- current+1
        ntemp <- as.integer(lines[current])
        ## print(ntemp)

        current <- current+1
        con <- textConnection(lines[current])
        temperatures <- as.double(read.delim(con, sep="", header=FALSE))
        close(con)

        con <- textConnection(lines[seq(ncoll)+current])
        rates <- read.delim(con, sep="", header=FALSE)
        close(con)
        rates <- rates[,-1]
        names(rates) <- c("up", "low", paste("rate.T", seq(ntemp), sep=""))

        collisions[[ipart]] <- list(partner=partner, T=temperatures, rates=rates)
        ## names(collisions)[ipart] <- partner
        
        current <- current+ncoll
    }
    
    moldata <- list(molecule=molecule, amass=amass,
                    levels=levels, transitions=transitions, collisions=collisions)
    moldata
}

