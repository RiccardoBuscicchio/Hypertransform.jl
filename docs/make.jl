using Documenter
using Hypertransform
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(
    sitename = "Hypertransform.jl",
    modules = [Hypertransform],
    pages = [
        "Home" => "index.md",
    ],
    plugins=[bib],
)

# Deploy the built docs (push to gh-pages)
deploydocs(
    repo = "github.com/RiccardoBuscicchio/Hypertransform.jl.git",
    branch = "gh-pages",
)
