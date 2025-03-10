inline float2 unity_voronoi_noise_randomVector4(float2 UV, float offset)
{
    float2x2 m = float2x2(15.27, 47.63, 99.41, 89.98);
    UV = frac(sin(mul(UV, m)) * 46839.32);
    return float2(sin(UV.y * offset) * 0.5 + 0.5, cos(UV.x * offset) * 0.5 + 0.5);
}
// 新しいボロノイノイズ生成関数
void Voronoi_Circle_float(float2 UV, float AngleOffset, float CellDensity, float CircleRadius, out float Out)
{
    float2 g = floor(UV * CellDensity);
    float2 f = frac(UV * CellDensity);
    float3 res = float3(8.0, 0.0, 0.0);

    // 元のボロノイノイズを使って各セルの位置を取得
    for(int y = -1; y <= 1; y++)
    {
        for(int x = -1; x <= 1; x++)
        {
            float2 lattice = float2(x, y);
            float2 offset = unity_voronoi_noise_randomVector4(lattice + g, AngleOffset);
            float2 basePos = lattice + offset;

            // 元のボロノイノイズの母点を除外し、黒丸の位置のみを使用する
            for(int i = 0; i < 9; i++)
            {
                float angle =  i * 2 * 3.14159 / 9.0; // 6等分
                float2 circlePos = basePos + float2(cos(angle), sin(angle)) * CircleRadius;
                float d = distance(circlePos, f);

                if(d < res.x)
                {
                    res = float3(d, circlePos.x, circlePos.y);
                }
            }
        }
    }
    Out = res.x;
}
