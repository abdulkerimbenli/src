"""
Find the routes according to the transit route solution file.
It extracts the route informatin given by the transit model

by July it is not working correctly

"""

using CSV, DataFrames, LightGraphs, SimpleWeightedGraphs, GraphPlot, DelimitedFiles

nNodes = 79
nLine = 8

# read the solutions file to get flows
df = CSV.read("/Users/akrm/Dropbox/Julia/transitComm/encapsulated/output/routes.csv", header = true, delim = ',',
                copycols = true)

# sort the solutions by commodity then by source
sort!(df, [:Route,:Source])

# commodityler üzerinden flowları toplamak için geçici dataframe
temp = DataFrame(Source = Int[], Sink = Int[], Route = Float64[], Flow = Float64[])

# ----------
# commodityler üzerinden akışları topla tüm edgelerdeki toplam akışı bul
# tüm networkdeki pathleri ayrı commodityler üzerinden değil
# bütün commodityler için bulur
for i in 1:size(df,1)
    floww = 0.0
    for j in 1:size(df,1)
        if df.Source[i] == df.Source[j] && df.Sink[i] == df.Sink[j] && df.Route[i] == df.Route[j]
            floww += df.flow[j]
        end
    end
    if floww != 0
        push!(temp, (df.Source[i], df.Sink[i], df.Route[i], floww))
    end
end
# aynı edgeleri sil
unique!(temp)


# ------- tek tek commodityler üzerinden pathlerin çıkarılması için kod parçası
# istenen commodity nin ayrıştırılması
for line = 1:nLine
    af = df
    af = af[af.Route .== line, :]

    # simple directed graph oluşturulması
    # 243 toplam node sayısını veriyor, datay göre değişir
    g = SimpleDiGraph(nNodes)

    # oluşturulan generic graph ın içine edge lerimizi ekliyoruz
    for i in 1:length(af.Source)
        add_edge!(g, af.Source[i], af.Sink[i])
    end

    # oluşan final graphın içerisindeki tüm simple pathler çıkarılır.
    # Ancak bu pathler tüm nodelar arasındak pathleri verir
    ListofPath = Any[]
    for path in enumerate_paths(dijkstra_shortest_paths(g, nNodes , allpaths=true))
        push!(ListofPath, path)
    end

    # kullanılmayan nodelardan dolayı boş pathler oluştuğu için bunlar sliniyor
    unique!(ListofPath)

    # tüm nodelardan tüm nodelara pathler çıkarılmıştı daha önceden
    # bize ise sadece source tan son sink e olan pathler lazım
    # ondan dolayı ara pathler silinir
    for i in 1:length(ListofPath), j in 1:length(ListofPath)
        if (i!= j) && issubset(ListofPath[i], ListofPath[j])
            ListofPath[i] = []
        end
    end
    unique!(ListofPath)

    # oluşan pathleri yazdır
    @show ListofPath

end
