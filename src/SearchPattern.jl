chain_generator(bond, targetAngle) = chain_generator(bond, targetAngle, targetAngle)

function chain_generator(bond, targetMin, targetMax, allowedMin = 0, allowedMax = 180)
    targetMin = deg2rad(targetMin)
    targetMax = deg2rad(targetMax)
    allowedMin = deg2rad(allowedMin)
    allowedMax = deg2rad(allowedMax)
    
    targetMean = (targetMin + targetMax) / 2
    targetMinCos = cos(targetMin)
    targetMaxCos = cos(targetMax)
    
    targetMeanSin = sin(targetMean)

    allowedMean = (allowedMin + allowedMax) / 2
    allowedMinCos = cos(allowedMin)
    allowedMaxCos = cos(allowedMax)
    
    targetArea = targetMinCos - targetMaxCos
    totalArea = allowedMinCos - allowedMaxCos

    relArea = targetArea / totalArea
    
    ΔϕMax = min(targetMean-allowedMin, allowedMax-targetMean)
    ΔϕFactor = totalArea/2targetMeanSin

    remainingSweep(tryFactor) = targetMean > allowedMean ?
        acos(allowedMaxCos + totalArea*tryFactor) :
        acos(allowedMinCos - totalArea*tryFactor)

    function(tryFactor)
        if tryFactor <= relArea
            ϕ = acos(targetMaxCos + targetArea*rand())
        elseif (k = tryFactor*ΔϕFactor) <= 1 && (Δϕ = asin(k)) <= ΔϕMax
            hemisphere = sign(rand() - sin(targetMean-Δϕ)/(2targetMeanSin*cos(Δϕ)))
            ϕ = targetMean + hemisphere * Δϕ
        else
            ϕ = remainingSweep(tryFactor)
        end
        bond, ϕ, 2pi*rand()
    end
end