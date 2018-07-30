using Documenter, UnitRoot

makedocs(modules=UnitRoot,
        doctest=true)

deploydocs(deps = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "https://github.com/p-chaim/UnitRoot.jl",
    julia  = "0.6.0",
    osname = "linux")
