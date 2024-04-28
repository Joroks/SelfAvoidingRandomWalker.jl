searchPattern(targetAngle) = searchPattern(targetAngle, targetAngle)

function searchPattern(targetMin, targetMax, allowedMin = 0, allowedMax = 180)
    targetMin = deg2rad(180 - targetMin)
    targetMax = deg2rad(180 - targetMax)
    allowedMin = deg2rad(180 - allowedMin)
    allowedMax = deg2rad(180 - allowedMax)
    
    targetMean = (targetMin + targetMax) / 2
    targetMinCos = cos(targetMin)
    targetMaxCos = cos(targetMax)
    
    targetMeanSin = sin(targetMean)

    allowedMean = (allowedMin + allowedMax) / 2
    allowedMinCos = cos(allowedMin)
    allowedMaxCos = cos(allowedMax)
    
    targetArea = targetMaxCos - targetMinCos
    totalArea = allowedMaxCos - allowedMinCos 

    relArea = targetArea / totalArea
    
    ΔϕMax = min(allowedMin-targetMean, targetMean-allowedMax)
    ΔϕFactor = totalArea/(sin(targetMin) + sin(targetMax))

    remainingSweep(tryFactor) = targetMean > allowedMean ?
        acos(allowedMinCos + totalArea*tryFactor) :
        acos(allowedMaxCos - totalArea*tryFactor)

    function randomDirection(tryFactor)
        if tryFactor <= relArea
            ϕ = acos(targetMinCos + targetArea*rand())
        elseif (k = tryFactor*ΔϕFactor) <= 1 &&  (Δϕ = asin(k)) <= ΔϕMax
            hemisphere = sign(rand() - sin(targetMean-Δϕ)/(2*targetMeanSin*cos(Δϕ)))
            ϕ = targetMean + hemisphere * Δϕ
        else
            ϕ = remainingSweep(tryFactor)
        end
        
        λ = 2*pi*rand()

        (sinλ, cosλ) = sincos(λ)
        (sinϕ, cosϕ) = sincos(ϕ)
    
        return SA[cosϕ 0 sinϕ; sinλ*sinϕ cosλ -cosϕ*sinλ; -cosλ*sinϕ sinλ cosλ*cosϕ]
    end

    return randomDirection
end