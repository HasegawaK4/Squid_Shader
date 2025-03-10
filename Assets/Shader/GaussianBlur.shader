Shader "Hidden/Custom/BlurPostProcess"
{
    Properties
    {
        // 必要に応じてプロパティを追加
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
            // カスタムポストプロセス処理を追加
            ENDCG
        }
    }
    Fallback "Hidden/Default"
}