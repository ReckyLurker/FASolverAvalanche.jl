using FASolverAvalanche
using Documenter

DocMeta.setdocmeta!(FASolverAvalanche, :DocTestSetup, :(using FASolverAvalanche); recursive=true)

makedocs(;
    modules=[FASolverAvalanche],
    authors="Tanish Jain",
    sitename="FASolverAvalanche.jl",
    format=Documenter.HTML(;
        canonical="https://ReckyLurker.github.io/FASolverAvalanche.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ReckyLurker/FASolverAvalanche.jl",
    devbranch="main",
)
