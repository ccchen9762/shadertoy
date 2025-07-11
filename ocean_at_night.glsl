const float PI = 3.1415926;

const int marchSteps = 64;
const int oceanOctaves = 4;
const int starOctaves = 4;

const mat2 rotation2D  = mat2(0.80, 0.60, -0.60, 0.80);
const mat2 rotation2Di = mat2(0.80, -0.60, 0.60, 0.80);

const mat2 rotation2D2  = mat2(2.4, 1.0, -1.0,  2.4);
const mat2 rotation2D2i = mat2(2.4, -1.0, 1.0,  2.4);

const vec3 oceanAmbeint = vec3(0.3, 0.6, 0.8) * 0.1;
const vec3 oceanDiffuse = vec3(0.3, 0.6, 0.8);
const vec3 oceanSpecular = vec3(1.0, 1.0, 1.0);
const float oceanHeight = 0.0;

const vec3 moonPosition = vec3(14.0, 10.0, -20.0);
float moonSize = 3.0;
const vec3 moonColor = vec3(1.0, 1.0, 0.85) * 0.6; 

const vec3 starColor = vec3(1.0, 1.0, 0.9) * 0.5;
const vec3 skyColor = vec3(0.35, 0.3, 0.5) * 0.5;

const vec3 lightDirection = normalize(vec3(0.5, 0.1, -1.0));
const vec3 lightColor = vec3(0.9, 0.9, 0.7) * 0.5;

const float minDistance = 0.01;
const float maxDistance = 50.0;


float hash(vec2 pos){
    float val = dot(pos, vec2(196.4, 548.9));
    return fract(sin(val) * 45125.29076);
}

// result .x = value noise .yz = x,z derivatives 
vec3 valueNoise(vec2 pos){
    vec2 i = floor(pos);
    vec2 f = fract(pos);
    vec2 u = f*f*(3.0-2.0*f); // 3x^2 - 2x^3
    vec2 du = 6.0*f*(1.0-f);  // 6x - 6x^2

    float a = hash(i + vec2(0.0, 0.0));
    float b = hash(i + vec2(1.0, 0.0));
    float c = hash(i + vec2(0.0, 1.0));
    float d = hash(i + vec2(1.0, 1.0));

    float k0 = a;
    float k1 = b-a;
    float k2 = c-a;
    float k3 = a-b-c+d;

    // bilinear interpolation (((d-c)*ux + c) - ((b-a)*ux + a))*uy + (b-a)*ux + a
    return vec3(-1.0 + 2.0 * (k0 + k1*u.x + k2*u.y + k3*u.x*u.y),
                2.0 * du * vec2(k1 + k3*u.y, k2 + k3*u.x));
}

float hash11(float p) {
    p = fract(p * .1031);
    p *= p + 33.33;
    p = fract(p * p);
    return p;
}

vec2 hash22(vec2 p) {
    vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.xx + p3.yz) * p3.zy);
}

float peakNoise(vec2 pos) {
    vec2 grid_scale_pos = pos * 7.0;

    vec2 i_pos = floor(grid_scale_pos);
    vec2 f_pos = fract(grid_scale_pos);

    float min_dist_to_peak = 1e5;
    float closest_peak_size_influence = 0.0; // Track size influence of the closest peak

    for (int y = -1; y <= 1; ++y) {
        for (int x = -1; x <= 1; ++x) {
            vec2 neighbor_cell = i_pos + vec2(float(x), float(y));
            vec2 offset_in_cell = hash22(neighbor_cell);
            float current_peak_size_factor = hash11(dot(neighbor_cell, vec2(5.3, 7.9)));

            vec2 peak_uv_absolute = (neighbor_cell + offset_in_cell) / 7.0;
            float dist_to_this_peak = length(pos - peak_uv_absolute);
            float scaled_dist = dist_to_this_peak * (0.8 + current_peak_size_factor * 0.7);

            if (scaled_dist < min_dist_to_peak) {
                min_dist_to_peak = scaled_dist;
            }
        }
    }

    float value = exp(-min_dist_to_peak * 80.0);
    return value;
}

vec3 oceanFBM(vec3 pos){
    float height = oceanHeight;
    float amplitude = 0.08;
    float gain = 0.3;
    float lacunarity = 1.85;
    vec2 derivative = vec2(0.0, 0.0);
    mat2 reverseMat = mat2(1.0, 0.0, 0.0, 1.0);

    vec2 uv = pos.xz;

    for(int i=0; i<oceanOctaves; i++){
        // add offset on each octaves to make the high frequency components move faster than lower freq
        vec2 animatedUV = uv + vec2(iTime * 0.6, iTime * 0.5);
        vec3 noise = valueNoise(animatedUV);
        height += amplitude*noise.x;                 // accumulate values		
        derivative += amplitude*reverseMat*noise.yz; // accumulate derivatives
        amplitude *= gain;
        uv = lacunarity*rotation2D*uv;
        reverseMat = lacunarity*rotation2Di*reverseMat;
    }

	return vec3(abs(pos.y - height), derivative);
}

float sdSphere(vec3 pos, float size){
  return length(pos)-size;
}

float diffuse(vec3 light, vec3 normal){
    return max(dot(light, normal), 0.0);
}
float specular(vec3 light, vec3 normal, vec3 rayDirection, float shiness){    
    return pow(max(dot(reflect(rayDirection, normal), light), 0.0), shiness);
}

vec3 getOceanColor(vec3 normal, vec3 rayDirection, float distance, float steepness){
    // foam effect
    vec3 foamOceanDiffuse = mix(oceanDiffuse, vec3(0.3, 0.4, 0.4), steepness * steepness * 40.0);
    
    vec3 diffuseColor = foamOceanDiffuse * diffuse(lightDirection, normal);
    vec3 specularColor = oceanSpecular * specular(lightDirection, normal, normalize(rayDirection), 50.0);

    float attenuation = 1.0 / (1.0 + distance * 0.05);
    vec3 color = oceanAmbeint + attenuation * (diffuseColor + specularColor * lightColor);
    
    return color;
}

vec3 getMoonColor(vec3 pos, vec3 normal, vec3 rayDirection) {
    float facingFactor = max(0.0, dot(normal, -rayDirection)); 
    float glowStrength = pow(facingFactor, 2.0);
    vec3 finalColor = moonColor * (0.5 + glowStrength * 0.5); 

    float detailNoise = 0.0;
    float scale = 10.0; 
    float persistence = 0.5;

    for (int i = 0; i < starOctaves; ++i) {
        detailNoise += peakNoise(pos.xz);
        scale *= 2.0;
        persistence *= 0.5;
    }

    detailNoise = detailNoise * 0.7 - 0.5;

    finalColor *= (1.0 + detailNoise * 0.2);

    finalColor = clamp(finalColor, 0.0, 1.5);

    return finalColor;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord){
    vec2 uv = (fragCoord.xy / iResolution.xy) * 2.0 - 1.0;
    uv.x *= iResolution.x / iResolution.y;
    
    // Camera setup
    vec3 cameraPosition  = vec3(0.0, 1.0, 3.0);
    vec3 cameraTarget = vec3(0.0, 1.75, 0.0);
    vec3 worldUp = vec3(0.0, 1.0, 0.0);
    
    vec3 cameraDirection = normalize(cameraTarget - cameraPosition);
    vec3 cameraRight = normalize(cross(cameraDirection, worldUp));
    vec3 cameraUp = cross(cameraRight, cameraDirection);
    
    float fov = radians(60.0);
    float focalLength = 1.0 / tan(fov * 0.5);
    
    vec3 rayOrigin = cameraPosition;
    vec3 rayDirection = normalize(
        cameraDirection * focalLength + 
        cameraRight * uv.x + 
        cameraUp * uv.y
    );
    
    vec3 pos = vec3(0.0);

    float totalDistance = 0.0; // distance travelled

    vec3 oceanMotion = vec3(0.0);
    float moonMotion = 0.0;
    float nextDistance = 0.0;
    
    vec3 color = vec3(0);

    for(int i=0; i<marchSteps; i++){
        pos = rayOrigin + rayDirection * totalDistance;

        oceanMotion = oceanFBM(pos);
        moonMotion = sdSphere(pos - moonPosition, moonSize);
        
        nextDistance = min(oceanMotion.x, moonMotion);
        
        if(nextDistance < minDistance)
            break;

        // map function is distance to surface at (x,z), not the actual shortest distance, need to optimize nextDistance
        totalDistance += nextDistance * 0.8 * (1.0 - 0.75 * min(length(oceanMotion.yz), 0.5));

        if(totalDistance > maxDistance){
            totalDistance = maxDistance + 1.0;
            break;
        }
    }
    
    //color = vec3(totalDistance * 0.05);
    if(totalDistance > maxDistance){      
        float noise = peakNoise(uv);
        float densityFactor = hash11(dot(uv, vec2(196.4, 548.9)));
        float starThreshold = 0.85 + densityFactor * -0.05;
        float starAppearance = step(starThreshold, noise);
        float time_factor = sin(hash(uv * iTime)) * 2.0 + 0.5;
        float brightFactor = hash(uv) + 0.1;
        vec3 finalStarColor = starColor * starAppearance * brightFactor * brightFactor * time_factor;

        color = vec3(skyColor - abs(rayDirection.y) * 0.3 + rayDirection.x * 0.05 + finalStarColor);
    }
    else{
        if(oceanMotion.x < moonMotion){
            // sea surface implicit function F(x,y,z)=y-h(x,z)=0
            // normal vector = ∇F = (-∂x, 1.0, -∂z);
            vec3 normal = normalize(vec3(-oceanMotion.y, 1.0, -oceanMotion.z));
            float steepness = length(vec2(oceanMotion.y, oceanMotion.z));
            color = getOceanColor(normal, rayDirection, totalDistance, steepness);
        }
        else{
            color = getMoonColor(pos, normalize(pos - moonPosition), rayDirection);
        }
    }

    fragColor = vec4(color, 1.0);
}