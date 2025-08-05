using Documenter
using Hypertransform

makedocs(
    sitename = "Hypertransform.jl",
    modules = [Hypertransform],
    pages = [
        "Home" => "index.md",
    ],
)

# Deploy the built docs (push to gh-pages)
deploydocs(
    repo = "github.com/RiccardoBuscicchio/Hypertransform.jl.git",
    branch = "gh-pages",
)
