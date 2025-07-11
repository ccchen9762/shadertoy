const int marchSteps = 64;
const int oceanOctaves = 5;

const mat2 rotation2D  = mat2(0.80, 0.60, -0.60, 0.80);
const mat2 rotation2Di = mat2(0.80, -0.60, 0.60, 0.80);

const mat2 rotation2D2  = mat2(1.6, 1.2, -1.2,  1.6);
const mat2 rotation2D2i = mat2(1.6, -1.2, 1.2,  1.6);

const vec3 oceanAmbeint = vec3(0.3, 0.6, 0.8) * 0.1;
const vec3 oceanDiffuse = vec3(0.3, 0.6, 0.8);
const vec3 oceanSpecular = vec3(1.0, 1.0, 1.0);
const float oceanHeight = 0.0;

const vec3 skyColor = vec3(0.35, 0.3, 0.5) * 0.5;

const vec3 lightDirection = normalize(vec3(0.35, 0.3, -0.6));
const vec3 lightColor = vec3(0.9, 0.8, 0.55) * 0.5;

const float minDistance = 0.001;
const float maxDistance = 30.0;


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

float diffuse(vec3 light, vec3 normal){
    return max(dot(light, normal), 0.0);
}
float specular(vec3 light, vec3 normal, vec3 rayDirection, float shiness){    
    return pow(max(dot(reflect(rayDirection, normal), light), 0.0), shiness);
}

vec3 getOceanColor(vec3 normal, vec3 rayDirection, float distance, float steepness){
    // foam effect
    vec3 foamOceanDiffuse = mix(oceanDiffuse, vec3(1.0, 1.0, 1.0), steepness * steepness * 50.0);
    
    vec3 diffuseColor = foamOceanDiffuse * diffuse(lightDirection, normal);
    vec3 specularColor = oceanSpecular * specular(lightDirection, normal, normalize(rayDirection), 10.0);
    
    float attenuation = 1.0 / (1.0 + distance * 0.2);
    vec3 color = oceanAmbeint + attenuation * (diffuseColor + specularColor * lightColor);
    
    
    return color;
}


void mainImage(out vec4 fragColor, in vec2 fragCoord){
    vec2 uv = (gl_FragCoord.xy / iResolution.xy) * 2.0 - 1.0;
    uv.x *= iResolution.x / iResolution.y;
    uv.y += 0.3; // viewing upward
    
    vec3 rayOrigin = vec3(0.0, 1.0, 3.0);
    vec3 rayDirection = normalize(vec3(uv, -1.0));
    vec3 pos = vec3(0.0);

    float totalDistance = 0.0; // distance travelled

    vec3 oceanMotion = vec3(0.0);
    float cloudMotion = 0.0;
    float nextDistance = 0.0;
    
    vec3 color = vec3(0);

    for(int i=0; i<marchSteps; i++){
        pos = rayOrigin + rayDirection * totalDistance;

        oceanMotion = oceanFBM(pos);
        
        nextDistance = oceanMotion.x;
        
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
        color = skyColor - abs(rayDirection.y) * 0.3;
    }
    else{
        // sea surface implicit function F(x,y,z)=y-h(x,z)=0
        // normal vector = ∇F = (-∂x, 1.0, -∂z);
        vec3 normal = normalize(vec3(-oceanMotion.y, 1.0, -oceanMotion.z));
        float steepness = length(vec2(oceanMotion.y, oceanMotion.z));
        color = getOceanColor(normal, rayDirection, totalDistance, steepness);
    }

    fragColor = vec4(color, 1.0);
}