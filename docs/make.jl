using Documenter
using Hypertransform

makedocs(
    sitename = "Hypertransform.jl",
    modules = [Hypertransform],
    pages = [
        "Home" => "index.md",
    ]
)
