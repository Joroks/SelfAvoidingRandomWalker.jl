module MakieExt

using SelfAvoidingRandomWalker
using GLMakie

function SelfAvoidingRandomWalker.visualize_search_pattern(n = 10000)
    fig = Figure()

    axes = GridLayout(fig[1, 1:2])

    ax1 = Axis3(axes[1,1];
        aspect = :data,
        limits = (-1..1, -1..1, -1..1)
    )

    ax2 = Axis(axes[1,2])

    target = IntervalSlider(fig[2,2], range=0:180); Label(fig[2,1], "target")
    allowed = IntervalSlider(fig[3,2], range=0:180); Label(fig[3,1], "allowed")
    visible = IntervalSlider(fig[4,2], range=1:n); Label(fig[4,1], "visible")

    k = range(0,1,n)
    points = lift(allowed.interval, target.interval) do allowed, target
        pattern = SelfAvoidingRandomWalker.searchPattern(target..., allowed...)
        map(1:n) do k
            try_factor = (k-rand())/n
            pattern(try_factor)[:,1]
        end |> stack
    end

    GLMakie.onany(points, visible.interval; update=true) do points, visible
        empty!(ax1); empty!(ax2)
        range = visible[1]:visible[2]
        scatter!(ax1, view(points, :, range), color=range)
        hist!(ax2, view(points, 1, range), color=:values, bins=100)
    end

    display(fig)
end

end