#==============================================================================
    Source code of Metrics and Heuristics

    Authors: 
    - Alex Paranahyba Abreu      <abreualexp@gmail.com>
    - Helio Yochihiro Fuchigami  <helio@dep.ufscar.br>
==============================================================================#


"""
Computation of the objective function.
"""
function IWT(p, sequencia)
    (n,m) = size(p)

    CT = zeros(Int, n, m)
    IT = zeros(Int, n, m)
    WT = zeros(Int, n, m)

    tarj = sequencia[1]
    CT[tarj,1] = p[tarj,1]

    for k=2:m
        CT[tarj,k] = CT[tarj,k-1] + p[tarj,k]
        WT[tarj,k-1] = 0
        IT[tarj,k] = sum(p[tarj,l] for l=1:k-1)
    end

    for j=2:n
        tari = sequencia[j-1]
        tarj = sequencia[j]
        CT[tarj,1] = CT[tari,1] + p[tarj,1]
    end

    for j=2:n
        for k=2:m
            tari = sequencia[j-1]
            tarj = sequencia[j]

            if ( CT[tarj,k-1] > CT[tari,k] )
                CT[tarj,k] = CT[tarj,k-1] + p[tarj,k]
                IT[tarj,k] = CT[tarj,k-1] - CT[tari,k]
            else
                CT[tarj,k] = CT[tari,k] + p[tarj,k]
                WT[tarj,k-1] = CT[tari,k] - CT[tarj,k-1]
            end
        end
    end
    IWT = sum(sum(IT[j,k] + WT[j,k] for j=1:n) for k=1:m)

    return IWT
end

#==============================================================================
    SCHEDULING METRICS
==============================================================================#

using Random
"""
Random sequence.
- Requires package 'Random'.
"""
function RAND(p)
    (n,m) = size(p)
    sequencia = Array{Int}(undef,n)

    for j = 1:n
        sequencia[j] = j
    end
    sequencia = shuffle!(sequencia)

    return sequencia
end

"""
Longest Processing Time (LPT) metric.
"""
function LPT(p)
    (n,m) = size(p)
    sequencia = Array{Int}(undef,n)
    ptotal = Array{Int}(undef,n)

    for j = 1:n
        ptotal[j] = sum(p[j,:])
        sequencia[j] = j
    end

    for i=1:n
        tari = sequencia[i]

        for j=i+1:n
            tarj = sequencia[j]

            if ( ptotal[tarj] > ptotal[tari] )
                sequencia[i] = tarj
                sequencia[j] = tari
                tari = tarj
            end
        end
    end
    return sequencia
end

"""
Shortest Processing Time (SPT) metric.
"""
function SPT(p)
    (n,m) = size(p)
    sequencia = Array{Int}(undef,n)
    ptotal = Array{Int}(undef,n)

    for j = 1:n
        ptotal[j] = sum(p[j,:])
        sequencia[j] = j
    end

    for i=1:n
        tari = sequencia[i]

        for j=i+1:n
            tarj = sequencia[j]

            if ( ptotal[tarj] < ptotal[tari] )
                sequencia[i] = tarj
                sequencia[j] = tari
                tari = tarj
            end
        end
    end

    return sequencia
end

"""
Shortest Front Idle Time (SFIT) metric.
"""
function SFIT(p)
    (n,m) = size(p)
    sequencia = Array{Int}(undef,n)
    ptotal = Array{Int}(undef,n)
    FI = Array{Int}(undef,n)

    for j = 1:n
        for x=1:m-1
            FI[j] = sum(p[j,k] for k=1:m-x)
        end
        sequencia[j] = j
    end

    for i=1:n
        tari = sequencia[i]

        for j=i+1:n
            tarj = sequencia[j]

            if ( FI[tarj] < FI[tari] )
                sequencia[i] = tarj
                sequencia[j] = tari
                tari = tarj
            end
        end
    end

    return sequencia
end

"""
Shortest Last Front Idle Time (SLFIT) metric.
"""
function SLFIT(p)
    (n,m) = size(p)
    sequencia = Array{Int}(undef,n)
    ptotal = Array{Int}(undef,n)
    FI = Array{Int}(undef,n)

    for j = 1:n
        FI[j] = sum(p[j,k] for k=1:m-1)
        sequencia[j] = j
    end

    for i=1:n
        tari = sequencia[i]

        for j=i+1:n
            tarj = sequencia[j]

            if ( FI[tarj] < FI[tari] )
                sequencia[i] = tarj
                sequencia[j] = tari
                tari = tarj
            end
        end
    end

    return sequencia
end

"""
Longest Front Idle Time (LFIT) metric.
"""
function LFIT(p)
    (n,m) = size(p)
    sequencia = Array{Int}(undef,n)
    ptotal = Array{Int}(undef,n)
    FI = Array{Int}(undef,n)

    for j = 1:n
        for x=1:m-1
            FI[j] = sum(p[j,k] for k=1:m-x)
        end
        sequencia[j] = j
    end

    for i=1:n
        tari = sequencia[i]

        for j=i+1:n
            tarj = sequencia[j]

            if ( FI[tarj] > FI[tari] )
                sequencia[i] = tarj
                sequencia[j] = tari
                tari = tarj
            end
        end
    end

    return sequencia
end

"""
Longest Last Front Idle Time (LLFIT) metric.
"""
function LLFIT(p)
    (n,m) = size(p)
    sequencia = Array{Int}(undef,n)
    ptotal = Array{Int}(undef,n)
    FI = Array{Int}(undef,n)

    for j = 1:n
        FI[j] = sum(p[j,k] for k=1:m-1)
        sequencia[j] = j
    end

    for i=1:n
        tari = sequencia[i]

        for j=i+1:n
            tarj = sequencia[j]

            if ( FI[tarj] > FI[tari] )
                sequencia[i] = tarj
                sequencia[j] = tari
                tari = tarj
            end
        end
    end

    return sequencia
end

"""
Shortest Back Idle Time (SBIT) metric.
"""
function SBIT(p)
    (n,m) = size(p)
    sequencia = Array{Int}(undef,n)
    ptotal = Array{Int}(undef,n)
    FI = Array{Int}(undef,n)

    for j = 1:n
        for x=2:m
            FI[j] = sum(p[j,k] for k=x:m)
        end
        sequencia[j] = j
    end

    for i=1:n
        tari = sequencia[i]

        for j=i+1:n
            tarj = sequencia[j]

            if ( FI[tarj] < FI[tari] )
                sequencia[i] = tarj
                sequencia[j] = tari
                tari = tarj
            end
        end
    end

    return sequencia
end

"""
Shortest First Back Idle Time (SFBIT) metric.
"""
function SFBIT(p)
    (n,m) = size(p)
    sequencia = Array{Int}(undef,n)
    ptotal = Array{Int}(undef,n)
    FI = Array{Int}(undef,n)

    # calcula tempo de processamento total de cada tarefa
    for j = 1:n
        FI[j] = sum(p[j,k] for k=2:m)
        sequencia[j] = j # inicializa o vetor sequencia[j]
    end

    # ordenação CRESCENTE pelo "Algoritmo de Seleção"
    for i=1:n
        tari = sequencia[i]

        # procura a tarefa com o maior tempo total de processamento após posição i
        for j=i+1:n
            tarj = sequencia[j]

            if ( FI[tarj] < FI[tari] )
                sequencia[i] = tarj
                sequencia[j] = tari
                tari = tarj
            end
        end
    end

    return sequencia
end

"""
Longest Front Idle Time (LBIT) metric.
"""
function LBIT(p)
    (n,m) = size(p)
    sequencia = Array{Int}(undef,n)
    ptotal = Array{Int}(undef,n)
    FI = Array{Int}(undef,n)

    for j = 1:n
        for x=2:m
            FI[j] = sum(p[j,k] for k=x:m)
        end
        sequencia[j] = j
    end

    for i=1:n
        tari = sequencia[i]

        for j=i+1:n
            tarj = sequencia[j]

            if ( FI[tarj] > FI[tari] )
                sequencia[i] = tarj
                sequencia[j] = tari
                tari = tarj
            end
        end
    end

    return sequencia
end

"""
Longest First Front Idle Time (LFBIT) metric.
"""
function LFBIT(p)
    (n,m) = size(p)
    sequencia = Array{Int}(undef,n)
    ptotal = Array{Int}(undef,n)
    FI = Array{Int}(undef,n)

    for j = 1:n
        FI[j] = sum(p[j,k] for k=2:m)
        sequencia[j] = j
    end

    for i=1:n
        tari = sequencia[i]

        for j=i+1:n
            tarj = sequencia[j]

            if ( FI[tarj] > FI[tari] )
                sequencia[i] = tarj
                sequencia[j] = tari
                tari = tarj
            end
        end
    end

    return sequencia
end


#==============================================================================
    HEURISTIC ALGORITHM
==============================================================================#

"""
NEH heuristic algorithm adapted from Nawaz, Enscore, and Ham (1983).
Implemented by HY Fuchigami and AP Abreu.
"""
function NEH_IW(p)
    (n, m) = size(p)
    sequencia = Array{Int}(undef,n)
    melhorseq = Array{Int}(undef,n)
    ptotal = Array{Int}(undef,n)
    CT = zeros(Int, n, m)
    IT = zeros(Int, n, m)
    WT = zeros(Int, n, m)


    # Gets initial sequence from the metric
    sequencia = RAND(p)
    #sequencia = LPT(p)
    #sequencia = SPT(p)
    #sequencia = SFIT(p)
    #sequencia = SLFIT(p)
    #sequencia = LFIT(p)
    #sequencia = LLFIT(p)
    #sequencia = SBIT(p)
    #sequencia = SFBIT(p)
    #sequencia = LBIT(p)
    #sequencia = LFBIT(p)

    melhorseq[1] = sequencia[1]

    for nparcial=2:n
        melhorIW = 1000000000000

        for j=1:nparcial-1
            sequencia[j] = melhorseq[j]
        end

        pos = nparcial
        while pos >= 1


            tarj = sequencia[1]
            CT[tarj,1] = p[tarj,1]

            for k=2:m
                CT[tarj,k] = CT[tarj,k-1] + p[tarj,k]
                WT[tarj,k-1] = 0
                IT[tarj,k] = sum(p[tarj,l] for l=1:k-1)
            end
            for j=2:nparcial
                tari = sequencia[j-1]
                tarj = sequencia[j]
                CT[tarj,1] = CT[tari,1] + p[tarj,1]
            end

            for j=2:nparcial
                for k=2:m
                    tari = sequencia[j-1]
                    tarj = sequencia[j]

                    if ( CT[tarj,k-1] > CT[tari,k] )
                        CT[tarj,k] = CT[tarj,k-1] + p[tarj,k]
                        IT[tarj,k] = CT[tarj,k-1] - CT[tari,k]
                    else
                        CT[tarj,k] = CT[tari,k] + p[tarj,k]
                        WT[tarj,k-1] = CT[tari,k] - CT[tarj,k-1]
                    end
                end
            end
            IWTparcial = sum(sum(IT[j,k] + WT[j,k] for j=1:n) for k=1:m)

            if ( IWTparcial < melhorIW )
                melhorIW = IWTparcial
                for j=1:nparcial
                    melhorseq[j] = sequencia[j]
                end
            end

            if ( pos > 1 )
                tarj = sequencia[pos]
                sequencia[pos] = sequencia[pos-1]
                sequencia[pos-1] = tarj
            end
            pos -= 1
        end
    end

    return melhorseq
end