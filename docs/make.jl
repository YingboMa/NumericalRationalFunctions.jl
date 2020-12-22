using NumericalRationalFunctions
using Documenter

makedocs(;
    modules=[NumericalRationalFunctions],
    authors="Yingbo Ma <mayingbo5@gmail.com> and contributors",
    repo="https://github.com/Yingbo Ma/NumericalRationalFunctions.jl/blob/{commit}{path}#L{line}",
    sitename="NumericalRationalFunctions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Yingbo Ma.github.io/NumericalRationalFunctions.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Yingbo Ma/NumericalRationalFunctions.jl",
)
