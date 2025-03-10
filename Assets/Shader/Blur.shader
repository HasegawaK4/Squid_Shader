Shader "Hidden/Custom/BlurShader"
{
    Properties
    {
        _BlurSize ("Blur Size", Float) = 1.0
    }

    SubShader
    {
        Tags { "RenderType"="Opaque" }
        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #include "Packages/com.unity.render-pipelines.high-definition/Runtime/RenderPipeline/ShaderLibrary/Core.hlsl"

            float _BlurSize;

            struct Attributes
            {
                float4 position : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct Varyings
            {
                float2 uv : TEXCOORD0;
                float4 position : SV_POSITION;
            };

            Varyings vert(Attributes v)
            {
                Varyings o;
                o.position = TransformObjectToHClip(v.position.xyz);
                o.uv = v.uv;
                return o;
            }

            half4 frag(Varyings i) : SV_Target
            {
                return half4(0, 0, 0, 1); // ここにブラー処理を追加
            }

            ENDCG
        }
    }

    Fallback "Hidden/Default"
}